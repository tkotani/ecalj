"""GPU locks shared by the ecalj jobs on one node: one lock file per GPU in /tmp/ecalj_res/.

A job takes the GPUs it will use with flock and hands only those to its programs in CUDA_VISIBLE_DEVICES.
flock is released by the kernel when the holder ends, so a job that crashes leaves no stale lock.  All or
nothing: a job that cannot get every GPU it wants holds none while it waits, so two jobs never hold one
GPU each and wait for each other.  The same directory is the place for later signals between jobs (none
yet).  2026-09-27.

    with GpuLock(want=2, label="gwsc hgw") as g:
        if g.devices is not None:
            env["CUDA_VISIBLE_DEVICES"] = g.visible

Environment:
    CUDA_VISIBLE_DEVICES   set by the user: only these GPUs are candidates (numbers only; otherwise no locking)
    ECALJ_GPU_LOCK=0       no locking at all
    ECALJ_GPU_WAIT=<sec>   give up after this many seconds (default: wait as long as it takes)
"""
from __future__ import annotations
import datetime
import fcntl
import os
import subprocess
import time

RES_DIR = "/tmp/ecalj_res"


def candidate_gpus() -> list[int] | None:
    """GPUs a job may use: CUDA_VISIBLE_DEVICES when set (None if it holds UUIDs or MIG names), else
    every GPU nvidia-smi lists ([] without GPU)."""
    v = os.environ.get("CUDA_VISIBLE_DEVICES")
    if v is not None and v.strip() != "":
        ids = [x.strip() for x in v.split(",") if x.strip() != ""]
        return [int(x) for x in ids] if all(x.isdigit() for x in ids) else None
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30).stdout
    except Exception:
        return []
    return list(range(sum(1 for line in out.splitlines() if line.startswith("GPU "))))


class GpuLock:
    """Hold `want` GPUs (at most as many as there are candidates) for the `with` block.
    devices is None when locking is off or not possible; the caller then leaves CUDA_VISIBLE_DEVICES alone."""

    def __init__(self, want: int = 1, label: str = "", wait: float | None = None):
        self.want = max(1, int(want or 1))
        self.label = label
        w = os.environ.get("ECALJ_GPU_WAIT")
        self.wait = wait if wait is not None else (float(w) if w else None)
        self.fds: list[int] = []
        self.devices: list[int] | None = None

    @property
    def visible(self) -> str:
        return ",".join(str(g) for g in self.devices or [])

    def __enter__(self):
        if os.environ.get("ECALJ_GPU_LOCK", "1") == "0":
            return self
        cands = candidate_gpus()
        if not cands:
            return self
        want = min(self.want, len(cands))
        os.makedirs(RES_DIR, exist_ok=True)
        try:
            os.chmod(RES_DIR, 0o1777)          # shared by every user of the node, like /tmp
        except OSError:
            pass
        t0 = time.time()
        told = 0.0
        while True:
            got = []
            for g in cands:
                fd = self._try(g)
                if fd is not None:
                    got.append((g, fd))
                    if len(got) == want:
                        break
            if len(got) == want:
                self.fds = [fd for _, fd in got]
                self.devices = [g for g, _ in got]
                self._note()
                return self
            for _, fd in got:                   # all or nothing
                os.close(fd)
            waited = time.time() - t0
            if self.wait is not None and waited >= self.wait:
                raise RuntimeError(f"GPU lock: {want} free GPU(s) among {cands} wanted; "
                                   f"in use: {self.holders(cands)}  (ECALJ_GPU_WAIT={self.wait:g}s)")
            if waited - told >= 600 or told == 0.0:
                print(f"GPU lock: waiting for {want} free GPU(s) among {cands}; in use: {self.holders(cands)}",
                      flush=True)
                told = max(waited, 1e-9)
            time.sleep(10 if self.wait is None else max(0.5, min(10.0, self.wait - waited)))

    def __exit__(self, *exc):
        for fd in self.fds:
            try:
                os.close(fd)                    # closing the file releases its flock
            except OSError:
                pass
        self.fds = []
        return False

    @staticmethod
    def _path(g: int) -> str:
        return f"{RES_DIR}/gpu{g}.lock"

    def _try(self, g: int) -> int | None:
        try:
            fd = os.open(self._path(g), os.O_RDWR | os.O_CREAT, 0o666)
        except PermissionError:
            fd = os.open(self._path(g), os.O_RDONLY)   # another user's file: flock works read-only too
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return fd
        except BlockingIOError:
            os.close(fd)
            return None

    def _note(self):
        """Who holds the GPU, for the message of a job that waits for it."""
        stamp = datetime.datetime.now().strftime("%m-%d %H:%M")
        for fd in self.fds:
            try:
                os.ftruncate(fd, 0)
                os.pwrite(fd, f"pid {os.getpid()} {self.label} since {stamp}\n".encode(), 0)
            except OSError:
                pass

    def holders(self, cands) -> str:
        out = []
        for g in cands:
            fd = self._try(g)
            if fd is not None:                  # free now; the file may still name an earlier holder
                os.close(fd)
                out.append(f"GPU {g}: free")
                continue
            try:
                with open(self._path(g)) as f:
                    txt = f.read().strip()
            except OSError:
                txt = ""
            out.append(f"GPU {g}: {txt or 'held'}")
        return "; ".join(out)
