from __future__ import annotations
from .utils import remove_files
import os
import shlex
import subprocess
import datetime
from pathlib import Path
from dataclasses import dataclass, field

import tomllib

START_TIME = datetime.datetime.now()
# Base directory of this module
CONFIG_DIR = Path(__file__).parent.parent

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
    """Build mpirun/srun command with simple placeholder replacement."""
    launcher = os.path.expandvars(cfg["launcher"])
    launcher_args_template = cfg.get("launcher_args", [])
    launcher_args = []

    for arg_template in launcher_args_template:
        s = arg_template
        if params.nprocs is not None:
            s = s.replace("{nprocs}", str(params.nprocs))
        if "{npernode}" in s:
            if params.npernode is None:
                continue
            s = s.replace("{npernode}", str(params.npernode))
        s = os.path.expandvars(s)
        launcher_args.extend(shlex.split(s))
    cmd_list = [launcher] + launcher_args

    if params.command:
        cmd_list.append(str(params.command))

    for arg in params.args:
        cmd_list.extend(shlex.split(os.path.expandvars(arg)))

    return cmd_list

def _build_env(cfg):
    """Build environment variables"""
    env = os.environ.copy()
    for k, v in cfg.get("env", {}).items():
        env[k] = v
    return env


# --- GPU slot coordination -------------------------------------------------
# gwsc runs *_gpu binaries (hvccfp0_mp_gpu, hgw_combined_mp_gpu, hsfp0_sc_mp_gpu,
# hrcxq_mp_gpu) directly via run_cmd. Without coordination, N parallel gwsc
# workers all pile onto GPU 0 -> CUDA_ERROR_OUT_OF_MEMORY (vcoulq_4 mkjp.f90:296).
# When a slot_scheduler_daemon is running we acquire a GPU slot here: it caps
# concurrency at the number of GPU slots and pins CUDA_VISIBLE_DEVICES to the
# slot index, so two GPU jobs land on GPU 0 and GPU 1 rather than colliding.
# No daemon -> no-op, so standalone use is unaffected.
_SLOT_SOCKET = "/tmp/slot_scheduler.sock"


class _GpuSlot:
    """Hold a GPU slot from slot_scheduler_daemon for the `with` block.
    slot_idx is None when no daemon is running (graceful fallback)."""

    def __init__(self, label="gwsc-gpu"):
        self.label = label
        self.sock = None
        self.slot_idx = None

    def __enter__(self):
        if not os.path.exists(_SLOT_SOCKET):
            return self
        try:
            import socket as _socket
            import json as _json
            s = _socket.socket(_socket.AF_UNIX, _socket.SOCK_STREAM)
            s.connect(_SLOT_SOCKET)
            fp = s.makefile("rwb")
            req = {"op": "REQUEST", "kind": "gpu",
                   "wid": self.label, "bin": self.label, "mpid": self.label}
            fp.write((_json.dumps(req) + "\n").encode())
            fp.flush()
            line = fp.readline()
            if not line:
                s.close()
                return self
            resp = _json.loads(line.decode())
            if resp.get("op") != "ASSIGN":
                s.close()
                return self
            self.sock = s
            self.slot_idx = resp.get("slot")
        except Exception:
            self.sock = None
            self.slot_idx = None
        return self

    def __exit__(self, *exc):
        if self.sock is not None:
            try:
                self.sock.close()  # daemon releases slot on disconnect
            except Exception:
                pass
            self.sock = None
        return False


def _needs_gpu_slot(command) -> bool:
    return command is not None and str(command).endswith("_gpu")


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
    gpu_ctx = _GpuSlot(label=f"runcmd:{Path(str(params.command)).name}") \
        if _needs_gpu_slot(params.command) else None
    if gpu_ctx is not None:
        gpu_ctx.__enter__()
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
            if gpu_ctx is not None and gpu_ctx.slot_idx is not None:
                env["CUDA_VISIBLE_DEVICES"] = str(gpu_ctx.slot_idx)
            dt = datetime.datetime.now() - START_TIME
            # Build the initial command for logging
            sec = dt.total_seconds()
            h = int(sec // 3600)
            m = int((sec % 3600) // 60)
            s = sec % 60
            stdout_info = f" stdout={stdout!r}" if stdout else ""
            stdin_info  = f" stdin={stdin_str!r}" if stdin_str else ""
            print(f"{h:02d}:{m:02d}:{s:06.3f}   {' '.join(map(str, cmd))}{stdout_info}{stdin_info}", end="", flush=True)
            t0 = datetime.datetime.now()
            # Execute the command
            rc = _run_mpi(cmd, env, stdin_str=stdin_str, stdout=out_stream)
            if rc == 0:
                elapsed_time = datetime.datetime.now() - t0
                sec = elapsed_time.total_seconds()
                if sec < 1:
                    print(f"  Elap. {sec:.3f}s", flush=True)
                elif sec < 60:
                    print(f"  Elap. {sec:.1f}s", flush=True)
                elif sec < 3600:
                    m = int(sec // 60)
                    s = int(sec % 60)
                    print(f"  Elap. {m}:{s:02d}", flush=True)
                else:
                    h = int(sec // 3600)
                    m = int((sec % 3600) // 60)
                    s = int(sec % 60)
                    print(f"  Elap. {h}:{m:02d}:{s:02d}", flush=True)
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
        if gpu_ctx is not None:
            gpu_ctx.__exit__(None, None, None)
        if out_stream:
            out_stream.close()
