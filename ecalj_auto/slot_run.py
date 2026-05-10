#!/usr/bin/env python3
"""
Run a command after acquiring a slot from the slot_scheduler daemon.

Usage:
  slot_run.py <kind> <wid> <bin_label> <mpid> -- <cmd> [args...]

Acquires CPU or GPU slot, runs the command (inherits stdin/stdout/stderr),
exits with the command's return code. Connection closes on exit → daemon
releases slot automatically.

For GPU acquisitions, sets CUDA_VISIBLE_DEVICES=<slot_idx> for the child.
"""
import json
import os
import socket
import subprocess
import sys

SOCKET_PATH = "/tmp/slot_scheduler.sock"


def main():
    args = sys.argv[1:]
    if "--" not in args or len(args) < 6:
        print("usage: slot_run.py <kind> <wid> <bin_label> <mpid> -- <cmd> [args...]",
              file=sys.stderr)
        sys.exit(2)
    sep = args.index("--")
    head = args[:sep]
    cmd = args[sep + 1:]
    if len(head) != 4:
        print("usage: slot_run.py <kind> <wid> <bin_label> <mpid> -- <cmd> [args...]",
              file=sys.stderr)
        sys.exit(2)
    kind, wid, bin_label, mpid = head

    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    sock.connect(SOCKET_PATH)
    f = sock.makefile("rwb")
    req = {"op": "REQUEST", "kind": kind, "wid": wid, "bin": bin_label, "mpid": mpid}
    f.write((json.dumps(req) + "\n").encode())
    f.flush()
    line = f.readline()
    if not line:
        print("daemon disconnected", file=sys.stderr)
        sys.exit(1)
    resp = json.loads(line.decode())
    if resp.get("op") != "ASSIGN":
        print(f"unexpected response: {resp}", file=sys.stderr)
        sys.exit(1)
    slot_idx = resp["slot"]

    env = os.environ.copy()
    if kind == "gpu":
        env["CUDA_VISIBLE_DEVICES"] = str(slot_idx)

    # Run command, inheriting stdin/stdout/stderr
    proc = subprocess.run(cmd, env=env)
    # On exit, sock is closed by Python runtime → daemon releases
    sock.close()
    sys.exit(proc.returncode)


if __name__ == "__main__":
    main()
