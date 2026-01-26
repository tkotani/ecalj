import os
import glob
import shutil
import shlex
import subprocess
import datetime
from pathlib import Path
from dataclasses import dataclass, field

import tomllib

START_TIME = datetime.datetime.now()
# Base directory of this module
CONFIG_DIR = Path(__file__).resolve().parent

# MPI execution parameters
@dataclass
class MPIParams:
    nprocs: int | None = 1 
    npernode: int | None = None
    command: str | Path | None = None
    args: list[str] = field(default_factory=list)


def _load_config(cluster_name):
    """Load cluster config"""
    path = CONFIG_DIR / "clusters.toml"
    if not path.is_file():
        raise FileNotFoundError(f"Cluster configuration file not found at: {path}")

    with open(path, "rb") as f:
        all_cfg = tomllib.load(f)
    
    if cluster_name in all_cfg:
        return all_cfg[cluster_name]

    # If not found, raise a more helpful error
    available = ", ".join(sorted(all_cfg.keys()))
    error_msg = (
        f"Cluster config not found: '{cluster_name}'.\n"
        f"Available clusters are: {available}"
    )
    raise RuntimeError(error_msg)


def _build_command(cfg, params: MPIParams) -> list[str]:
    """Build mpirun/srun command"""
    launcher = cfg["launcher"]
    launcher_args_template = cfg.get("launcher_args", [])
    launcher_args = []
    format_kwargs = {"nprocs": params.nprocs}
    if params.npernode is not None:
        format_kwargs["npernode"] = params.npernode

    for arg_template in launcher_args_template:
        # if npernode is required but not provided, skip this template
        if "{npernode}" in arg_template and params.npernode is None:
            continue

        formatted = arg_template.format(**format_kwargs)
        launcher_args.extend(shlex.split(formatted))
    return [launcher] + launcher_args + [params.command] + params.args


def _build_env(cfg):
    """Build environment variables"""
    env = os.environ.copy()
    for k, v in cfg.get("env", {}).items():
        env[k] = v
    return env


def _run_mpi(cmd, env, stdin_str=None, stdout=None):
    """Execute MPI command"""
    kwargs = dict(env=env, stdout=stdout, text=True)
    if stdin_str is not None:
        kwargs["input"] = stdin_str
    proc = subprocess.run(cmd, **kwargs)
    return proc.returncode


def run_cmd(cluster: str,
            params: MPIParams,
            retry: bool = False,
            stdin_str: str | None = None,
            stdout: str | None = None) -> None:
    """Execute Command with Retry Logic"""
    cluster = cluster or "default"
    cfg = _load_config(cluster)
    out_stream = open(stdout, "w") if stdout else None
    try:
        n = params.nprocs
        pnode = params.npernode
        while True:
            # Build MPIParams for this attempt
            p = MPIParams(
                nprocs=n,
                npernode=pnode,
                command=params.command,
                args=params.args,
            )
            # Build command and environment for this attempt
            cmd = _build_command(cfg, p)
            cmd = [str(x) for x in cmd]
            env = _build_env(cfg)
            dt = datetime.datetime.now() - START_TIME
            # Build the initial command for logging
            redir = f" > {stdout}" if stdout else ""
            print(f"{dt}   {' '.join(map(str, cmd))}{redir}", flush=True)
            # Execute the command
            rc = _run_mpi(cmd, env, stdin_str=stdin_str, stdout=out_stream)
            if rc == 0:
                return  # Success
            # If retry is disabled, fail immediately
            if not retry:
                raise RuntimeError(f"Command failed: {' '.join(map(str, cmd))}")
            # Retry mode: shrink nprocs and npernode
            if n == 1:
                raise RuntimeError("MPI failed even with nprocs=1")
            n = max(1, n // 2)
            if pnode is not None:
                pnode = max(1, pnode // 2)
            print(f"Retrying with nprocs={n}, npernode={pnode}")
    finally:
        if out_stream:
            out_stream.close()


def _read_bmix_from_ctrl(target: str) -> float:
    """Reads BMIX/b value from ctrl.{target} file."""
    import re
    ctrl_file = f'ctrl.{target}'
    if not Path(ctrl_file).is_file():
        raise FileNotFoundError(f"Control file not found: {ctrl_file}")
    with open(ctrl_file, 'r') as f:
        text = f.read()
    # Search for BMIX or b value
    bval_search = re.search(r'BMIX\s*=\s*([0-9.]+)', text, re.I)
    if not bval_search:
        bval_search = re.search(r'\bb\s*=\s*([0-9.]+)', text, re.I)
    if not bval_search:
        raise ValueError(f"Cannot find BMIX or b value in {ctrl_file}")

    return float(bval_search.group(1))


def _check_lmf_convergence(save_file: str) -> tuple[bool, str]:
    """
    Checks the save file for convergence markers.
    Returns a tuple of (converged_status, last_line).
    """
    with open(save_file, 'r') as f:
        lines = f.readlines()
    last_line = lines[-1].strip()
    if not last_line and len(lines) > 1:
        last_line = lines[-2].strip()
    if last_line:
        first_word = last_line.split()[0]
        if first_word.lower() in ['c', 'x', 'done']:
            return True, last_line
    return False, last_line


def _prepare_for_lmf_retry(rst_file: str):
    """Cleans up files and restores rst file for a retry."""
    # Clean up mix files
    for file in glob.glob("mix*") + glob.glob("__mix*"):
        Path(file).unlink()
    # Restore rst file from backup
    rst_bk_file = f'{rst_file}.bk'
    if Path(rst_bk_file).exists():
        if Path(rst_file).exists():
            Path(rst_file).unlink()
        shutil.move(rst_bk_file, rst_file)


const_b = {}
def run_lmf(cluster: str,
            target: str,
            params: MPIParams,
            bmix_reduction: bool = False,
            stdin_str: str | None = None,
            stdout: str | None = None) -> None:
    """
    Executes the 'lmf' command with special handling for convergence.
    Wraps `run_cmd` to automatically reduce `bmix` on convergence failure.
    """
    bval = _read_bmix_from_ctrl(target)
    rst_file = f'rst.{target}'
    save_file = f'save.{target}'

    if params.command in const_b:
        bval = const_b[params.command]
        print(f"Using cached b-value {bval} for {params.command}", flush=True)

    while True:
        if Path(rst_file).is_file():
            shutil.copy(rst_file, f'{rst_file}.bk')
        current_params = MPIParams(
            nprocs=params.nprocs,
            npernode=params.npernode,
            command=params.command,
            args=params.args + [f'-vb={bval}']
        )
        try:
            run_cmd(cluster, current_params, retry=False, stdin_str=stdin_str, stdout=stdout)
            converged, last_line = _check_lmf_convergence(save_file)
            if converged:
                print(f"lmf successful with b={bval}", flush=True)
                const_b[params.command] = bval
                if Path(f'{rst_file}.bk').exists():
                    Path(f'{rst_file}.bk').unlink()
                return
            raise RuntimeError(f'lmf did not converge. Last line of {save_file}: "{last_line}"')
        except RuntimeError as e:
            print(f'lmf failed with b={bval}: {e}', flush=True)
            if not bmix_reduction:
                raise RuntimeError(f"lmf failed with b={bval} and bmix_reduction is off.") from e
            _prepare_for_lmf_retry(rst_file)
            bval = round(bval - 0.05, 2)
            if bval < 0.05:
                raise RuntimeError("lmf failed even after reducing bmix to minimum.") from e
            print(f'Retrying with b={bval}', flush=True)
