#!/usr/bin/env python3
"""GW1500 slot scheduler daemon (FIFO).

Manages 2 CPU slots + 2 GPU slots via unix socket /tmp/slot_scheduler.sock.

Protocol (JSON-line, one message per line):
  Client -> Daemon:
    {"op":"REQUEST","kind":"cpu"|"gpu","wid":...,"bin":...,"mpid":...}
    {"op":"STATUS"}
  Daemon -> Client:
    {"op":"ASSIGN","slot":0|1}
    {"op":"STATUS","slots":{...},"queues":{...}}
"""
import asyncio
import json
import os
import sys
import time
from collections import deque

SOCKET_PATH = "/tmp/slot_scheduler.sock"

N_SLOTS = {"cpu": 2, "gpu": 2}
slots = {(kind, i): None for kind in N_SLOTS for i in range(N_SLOTS[kind])}
queues = {"cpu": deque(), "gpu": deque()}


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"{ts} {msg}\n")
    sys.stderr.flush()


def find_free_slot(kind):
    for i in range(N_SLOTS[kind]):
        if slots[(kind, i)] is None:
            return i
    return None


def assign(kind, idx, req, writer):
    slots[(kind, idx)] = (req["wid"], req.get("mpid", "?"), req.get("bin", "?"), writer)
    wid = req["wid"]
    mpid = req.get("mpid", "?")
    binname = req.get("bin", "?")
    log(f"ASSIGN {kind}_slot_{idx} -> {wid} {mpid} {binname}")


async def send_assign(writer, idx):
    try:
        msg = json.dumps({"op": "ASSIGN", "slot": idx}) + "\n"
        writer.write(msg.encode())
        await writer.drain()
    except Exception as e:
        log(f"send_assign error: {e}")


def release_by_writer(writer):
    for key, holder in slots.items():
        if holder is not None and holder[3] is writer:
            wid, mpid, binname, _ = holder
            slots[key] = None
            kind, idx = key
            log(f"RELEASE {kind}_slot_{idx} <- {wid} {mpid} {binname}")
            return kind, idx
    return None


async def handle_client(reader, writer):
    try:
        line = await reader.readline()
        if not line:
            return
        req = json.loads(line.decode().strip())
        op = req.get("op")

        if op == "STATUS":
            slots_snap = {}
            for (k, i), h in slots.items():
                slots_snap[f"{k}_{i}"] = (h[0], h[1], h[2]) if h else None
            queues_snap = {}
            for kind, q in queues.items():
                queues_snap[kind] = [(r["wid"], r.get("mpid", "?"), r.get("bin", "?")) for r, _ in q]
            writer.write((json.dumps({"op": "STATUS", "slots": slots_snap, "queues": queues_snap}) + "\n").encode())
            await writer.drain()
            return

        if op != "REQUEST":
            log(f"unknown op: {op}")
            return

        kind = req.get("kind")
        if kind not in N_SLOTS:
            log(f"bad kind: {kind}")
            return

        idx = find_free_slot(kind)
        if idx is not None:
            assign(kind, idx, req, writer)
            await send_assign(writer, idx)
        else:
            wid = req["wid"]
            mpid = req.get("mpid", "?")
            binname = req.get("bin", "?")
            log(f"QUEUE {kind} <- {wid} {mpid} {binname} (queue size {len(queues[kind])+1})")
            queues[kind].append((req, writer))

        # Wait until client disconnects (= release signal)
        while True:
            data = await reader.read(1024)
            if not data:
                break
    except Exception as e:
        log(f"handler error: {e}")
    finally:
        # Cleanup: remove from queue if waiting
        for kind, q in queues.items():
            for entry in list(q):
                if entry[1] is writer:
                    q.remove(entry)
                    log(f"DEQUEUE {kind} <- {entry[0]['wid']} (disconnect before assign)")
        # Release any held slot
        rel = release_by_writer(writer)
        if rel is not None:
            kind, _ = rel
            if queues[kind]:
                next_req, next_writer = queues[kind].popleft()
                free_idx = find_free_slot(kind)
                if free_idx is not None:
                    assign(kind, free_idx, next_req, next_writer)
                    await send_assign(next_writer, free_idx)
        try:
            writer.close()
            await writer.wait_closed()
        except Exception:
            pass


async def main():
    try:
        os.unlink(SOCKET_PATH)
    except FileNotFoundError:
        pass
    server = await asyncio.start_unix_server(handle_client, path=SOCKET_PATH)
    os.chmod(SOCKET_PATH, 0o666)
    log(f"slot_scheduler_daemon listening on {SOCKET_PATH}, PID={os.getpid()}")
    async with server:
        await server.serve_forever()


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        log("interrupted")
