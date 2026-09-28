#!/bin/bash
# Keep the ecalj source identical on every machine we compute on.
#
#   sync_ecalj_src.sh <host> [<remote dir>]      ship this HEAD to <host> (default dir: see default_dir below;
#                                                on kt1 the development tree ~/ecalj_dev, never the production ~/ecalj)
#   sync_ecalj_src.sh --samples <host> [...]     the whole tracked tree, Samples too (testecalj --all runs in
#                                                Samples/TestInstall; a new machine, 2026-09-28 kr7)
#   ALLOW_DIRTY=1 sync_ecalj_src.sh ...          ship HEAD although SRC has uncommitted changes (they are NOT shipped)
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
#   1) the archive includes SRC/CMakeLists.txt, so cmake re-configures.  REBUILD WITH THE
#      PROJECT INSTALLER, not a bare make:
#        cd ~/ecalj_dev && python3 InstallAll.py --fc nvfortran --gpu --gemmul8 --bindir ~/bin_dev   (kt1)
#      SRC/CMakeLists.txt selects the flag set by MATCHING THE STRING IN $FC, and
#      BUILD_DIR is SRC/build_$FC, so FC must be the bare compiler name.  A bare
#      `make` without FC stops at "Fortran compiler must be set via FC" and leaves the
#      old .so in place; worse, FC=<path>/mpifort silently matches the "ifort" branch
#      (m-p-*i-f-o-r-t*) and hands nvfortran Intel flags (-assume, -init:snan).  Both
#      failures leave a stale library that looks like a successful sync.
#      InstallAll.py also caps parallelism at 8 (nvfortran ICEs above that).
#   2) confirm the change actually reached the library, e.g.
#        strings SRC/build_nvfortran/libecaljF.so | grep -c ZmloSig
#      The executables are thin wrappers; the code lives in libecaljF*.so.
#
# It deliberately does NOT use git on the remote side.  kt1 and kr5 are separate
# checkouts (kt1 even carries a commit local does not), so a pull/merge there is a
# conflict waiting to happen.  The local repo is the single source of truth; a remote
# whose .ecalj_rev differs is simply stale and gets overwritten.
set -u
HOSTS_DEFAULT="kt1 kr5 kr7"
default_dir(){  # 2026-09-28: kt1's /home/takao/ecalj is the production tree (-> ~/bin) and must not be overwritten by a sync
  case $1 in kt1) echo /home/takao/ecalj_dev ;; *) echo /home/takao/ecalj ;; esac
}
LOCAL_DIR=$(cd "$(dirname "$0")/.." && pwd)

rev_local(){ (cd "$LOCAL_DIR" && git rev-parse --short HEAD); }
dirty_local(){ (cd "$LOCAL_DIR" && git status --porcelain -- SRC InstallAll.py | head -5); }

check_host(){  # $1=host $2=dir
  local h=$1 d=${2:-$(default_dir $1)}
  local r
  r=$(timeout 25 ssh -o ConnectTimeout=8 -o BatchMode=yes "$h" "cat $d/SRC/.ecalj_rev 2>/dev/null | head -1" 2>/dev/null)
  if [ -z "$r" ]; then
    r=$(timeout 25 ssh -o ConnectTimeout=8 -o BatchMode=yes "$h" "echo up" 2>/dev/null)
    [ -z "$r" ] && { printf '%-6s %s\n' "$h" "UNREACHABLE"; return; }
    printf '%-6s %s\n' "$h" "no .ecalj_rev (never synced by this script)"
    return
  fi
  printf '%-6s %s\n' "$h" "$r  ($d)"
}

ship(){  # $1=host $2=dir
  local h=$1 d=${2:-$(default_dir $1)}
  local rev tmp tarf
  rev=$(rev_local)
  local dirt; dirt=$(dirty_local)
  if [ -n "$dirt" ] && [ "${ALLOW_DIRTY:-0}" != 1 ]; then
    echo "REFUSING: SRC has uncommitted changes -- commit first so the marker means something:"
    echo "$dirt"; return 1
  fi
  [ -n "$dirt" ] && echo "NOTE: SRC has uncommitted changes; they are NOT shipped (HEAD $rev only)"
  local paths="SRC InstallAll.py"; [ "$SAMPLES" = 1 ] && paths=""      # "": the whole tracked tree
  tmp=$(mktemp -d); tarf=$tmp/ecalj_src.tar
  ( cd "$LOCAL_DIR" && git archive --format=tar HEAD $paths > "$tarf" ) || return 1
  mkdir -p "$tmp/SRC"
  printf '%s  %s  from %s\n' "$rev" "$(date '+%Y-%m-%d %H:%M')" "$(hostname)" > "$tmp/SRC/.ecalj_rev"
  ( cd "$tmp" && tar rf "$tarf" SRC/.ecalj_rev ) || return 1
  gzip -1 -f "$tarf"
  echo "shipping $rev to $h:$d  ($(du -h "$tarf.gz" | cut -f1))"
  timeout 900 scp -q "$tarf.gz" "$h:/tmp/ecalj_src.tar.gz" || { echo "scp FAILED"; rm -rf "$tmp"; return 1; }
  timeout 300 ssh -o ConnectTimeout=10 "$h" "mkdir -p $d && tar xzf /tmp/ecalj_src.tar.gz -C $d && rm -f /tmp/ecalj_src.tar.gz && cat $d/SRC/.ecalj_rev"
  rm -rf "$tmp"
}

SAMPLES=0; [ "${1:-}" = --samples ] && { SAMPLES=1; shift; }
case "${1:-}" in
  --check-all)
    printf '%-6s %s\n' local "$(rev_local)  (HEAD)"
    d=$(dirty_local); [ -n "$d" ] && echo "       local SRC is DIRTY:" && echo "$d"
    for h in $HOSTS_DEFAULT; do check_host "$h"; done ;;
  --check) shift; check_host "$@" ;;
  "") echo "usage: $0 <host> [<remote dir>] | --check <host> | --check-all"; exit 1 ;;
  *) ship "$@" ;;
esac
