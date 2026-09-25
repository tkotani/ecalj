#!/bin/bash
# Keep the ecalj source identical on every machine we compute on.
#
#   sync_ecalj_src.sh <host> [<remote dir>]      ship this HEAD to <host>
#   sync_ecalj_src.sh --check <host> [...]       report what <host> has
#   sync_ecalj_src.sh --check-all                report local + every known host
#
# Ships the TRACKED tree of the current HEAD (SRC + InstallAll.py) plus a marker
# file SRC/.ecalj_rev holding the hash and the date, so each machine can say which
# revision its binaries were built from.  Written 2026-09-25 after a whole morning of
# LiTi2O4 runs turned out to have been made with binaries that predated the fix under
# test: kt1's m_sigmlo.f90 had none of it, and nothing in the workflow would have said so.
#
# AFTER SYNCING you must rebuild, and the marker does NOT prove the binaries changed:
#   1) the archive includes SRC/CMakeLists.txt, so cmake re-configures, and this project
#      refuses to configure without FC.  On kt1:
#        cd ~/ecalj/SRC/build_nvfortran && \
#        FC=/opt/nvidia/hpc_sdk/Linux_x86_64/2026/comm_libs/mpi/bin/mpifort \
#        make -j16 lmf lmf_mp_gpu mlo hqpe_sc
#      Without FC it stops at "Fortran compiler must be set via FC" and leaves the OLD
#      .so in place -- a silent no-op that looks like a successful sync.
#   2) confirm the change actually reached the library, e.g.
#        strings SRC/build_nvfortran/libecaljF.so | grep -c ZmloSig
#      The executables are thin wrappers; the code lives in libecaljF*.so.
#
# It deliberately does NOT use git on the remote side.  kt1 and kr5 are separate
# checkouts (kt1 even carries a commit local does not), so a pull/merge there is a
# conflict waiting to happen.  The local repo is the single source of truth; a remote
# whose .ecalj_rev differs is simply stale and gets overwritten.
set -u
HOSTS_DEFAULT="kt1 kr5"
REMOTE_DIR_DEFAULT=/home/takao/ecalj
LOCAL_DIR=$(cd "$(dirname "$0")/.." && pwd)

rev_local(){ (cd "$LOCAL_DIR" && git rev-parse --short HEAD); }
dirty_local(){ (cd "$LOCAL_DIR" && git status --porcelain -- SRC InstallAll.py | head -5); }

check_host(){  # $1=host $2=dir
  local h=$1 d=${2:-$REMOTE_DIR_DEFAULT}
  local r
  r=$(timeout 25 ssh -o ConnectTimeout=8 -o BatchMode=yes "$h" "cat $d/SRC/.ecalj_rev 2>/dev/null | head -1" 2>/dev/null)
  if [ -z "$r" ]; then
    r=$(timeout 25 ssh -o ConnectTimeout=8 -o BatchMode=yes "$h" "echo up" 2>/dev/null)
    [ -z "$r" ] && { printf '%-6s %s\n' "$h" "UNREACHABLE"; return; }
    printf '%-6s %s\n' "$h" "no .ecalj_rev (never synced by this script)"
    return
  fi
  printf '%-6s %s\n' "$h" "$r"
}

ship(){  # $1=host $2=dir
  local h=$1 d=${2:-$REMOTE_DIR_DEFAULT}
  local rev tmp tarf
  rev=$(rev_local)
  local dirt; dirt=$(dirty_local)
  if [ -n "$dirt" ]; then
    echo "REFUSING: SRC has uncommitted changes -- commit first so the marker means something:"
    echo "$dirt"; return 1
  fi
  tmp=$(mktemp -d); tarf=$tmp/ecalj_src.tar
  ( cd "$LOCAL_DIR" && git archive --format=tar HEAD SRC InstallAll.py > "$tarf" ) || return 1
  mkdir -p "$tmp/SRC"
  printf '%s  %s  from %s\n' "$rev" "$(date '+%Y-%m-%d %H:%M')" "$(hostname)" > "$tmp/SRC/.ecalj_rev"
  ( cd "$tmp" && tar rf "$tarf" SRC/.ecalj_rev ) || return 1
  gzip -1 -f "$tarf"
  echo "shipping $rev to $h:$d  ($(du -h "$tarf.gz" | cut -f1))"
  timeout 900 scp -q "$tarf.gz" "$h:/tmp/ecalj_src.tar.gz" || { echo "scp FAILED"; rm -rf "$tmp"; return 1; }
  timeout 300 ssh -o ConnectTimeout=10 "$h" "mkdir -p $d && tar xzf /tmp/ecalj_src.tar.gz -C $d && rm -f /tmp/ecalj_src.tar.gz && cat $d/SRC/.ecalj_rev"
  rm -rf "$tmp"
}

case "${1:-}" in
  --check-all)
    printf '%-6s %s\n' local "$(rev_local)  (HEAD)"
    d=$(dirty_local); [ -n "$d" ] && echo "       local SRC is DIRTY:" && echo "$d"
    for h in $HOSTS_DEFAULT; do check_host "$h"; done ;;
  --check) shift; check_host "$@" ;;
  "") echo "usage: $0 <host> [<remote dir>] | --check <host> | --check-all"; exit 1 ;;
  *) ship "$@" ;;
esac
