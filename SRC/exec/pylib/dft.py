import re
import glob
import shutil
from pathlib import Path

from .utils import remove_files
from .run_cmd import run_cmd, MPIParams


def _read_bmix_from_ctrl(target: str) -> float:
    toml_file = f'ctrlg.{target}.toml'
    if Path(toml_file).is_file():
        import tomllib
        with open(toml_file, 'rb') as f:
            data = tomllib.load(f)
        b = data.get('iter', {}).get('b')
        if b is not None:
            return float(b)
        raise ValueError(f"Cannot find [iter] b in {toml_file}")
    ctrl_file = f'ctrl.{target}'
    if not Path(ctrl_file).is_file():
        raise FileNotFoundError(f"Control file not found: neither {toml_file} nor {ctrl_file}")
    with open(ctrl_file, 'r') as f:
        text = f.read()
    bval_search = re.search(r'BMIX\s*=\s*([0-9.]+)', text, re.I)
    if not bval_search:
        bval_search = re.search(r'\bb\s*=\s*([0-9.]+)', text, re.I)
    if not bval_search:
        raise ValueError(f"Cannot find BMIX or b value in {ctrl_file}")
    return float(bval_search.group(1))


def _check_lmf_convergence(save_file: str) -> tuple[bool, str]:
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
    for file in glob.glob("mix*") + glob.glob("__mix*"):
        Path(file).unlink()
    rst_bk_file = f'{rst_file}.bk'
    if Path(rst_bk_file).exists():
        if Path(rst_file).exists():
            Path(rst_file).unlink()
        shutil.move(rst_bk_file, rst_file)


def _ensure_ctrl(target):
    if not Path(f"ctrlg.{target}.toml").is_file() and not Path(f"ctrl.{target}").is_file():
        raise RuntimeError(f"No ctrl file (neither ctrlg.{target}.toml nor ctrl.{target})")


_const_b: dict = {}


def run_lmf(cluster: str,
            target: str,
            params: MPIParams,
            bmix_reduction: bool = False,
            stdin_str: str | None = None,
            stdout: str | None = None) -> None:
    """Run lmf with automatic bmix reduction on convergence failure."""
    _ensure_ctrl(target)
    bval = _read_bmix_from_ctrl(target)
    rst_file = f'rst.{target}'
    save_file = f'save.{target}'

    if params.command in _const_b:
        bval = _const_b[params.command]
        print(f"Using cached b-value {bval} for {params.command}", flush=True)

    while True:
        if Path(rst_file).is_file():
            shutil.copy(rst_file, f'{rst_file}.bk')
        current_params = MPIParams(
            nprocs=params.nprocs,
            npernode=params.npernode,
            command=params.command,
            args=params.args + [f'--ctrlg:iter.b={bval}']
        )
        try:
            run_cmd(cluster, current_params, retry=False, stdin_str=stdin_str, stdout=stdout)
            converged, last_line = _check_lmf_convergence(save_file)
            if converged:
                print(f"lmf successful with b={bval}", flush=True)
                _const_b[params.command] = bval
                if Path(f'{rst_file}.bk').exists():
                    Path(f'{rst_file}.bk').unlink()
                return
            raise RuntimeError(f'lmf did not converge. Last line of {save_file}: "{last_line}"')
        except RuntimeError as e:
            print(f'lmf failed with b={bval}: {e}', flush=True)
            if not bmix_reduction:
                raise RuntimeError(f"lmf failed with b={bval} and bmix_reduction is off.") from e
            _prepare_for_lmf_retry(rst_file)
            remove_files("__mixm")
            bval = round(bval - 0.05, 2)
            if bval < 0.05:
                raise RuntimeError("lmf failed even after reducing bmix to minimum.") from e
            print(f'Retrying with b={bval}', flush=True)


def cal_dft(cluster: str, target: str, mpi_size: int, exec_dir: Path, extra_args: list | None = None) -> None:
    """Run standard DFT: lmfa (atomic SCF) then lmf self-consistent field."""
    extra = extra_args or []
    run_cmd(cluster, MPIParams(nprocs=1, command=exec_dir / "lmfa",
                               args=[target, *extra]), stdout="llmfa")
    run_lmf(cluster, target, MPIParams(nprocs=mpi_size, command=exec_dir / "lmf",
                                       args=[target, *extra]),
            bmix_reduction=True, stdout="llmf")
