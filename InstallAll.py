#!/usr/bin/env python3
import os
import sys

if sys.version_info < (3, 11):
    sys.stderr.write(
        "ERROR: ecalj requires Python 3.11 or newer (found {}.{}.{}).\n"
        "       The build helper and the testecalj scripts use stdlib `tomllib`\n"
        "       and `contextlib.chdir`, both added in Python 3.11.\n"
        "\n"
        "       Install a recent Python locally and re-run this script. Examples:\n"
        "         pyenv install 3.12.13 && pyenv global 3.12.13\n"
        "         curl -LsSf https://astral.sh/uv/install.sh | sh \\\n"
        "             && uv python install 3.12 \\\n"
        "             && ln -sf \"$(uv python find 3.12)\" ~/.local/bin/python3\n"
        "\n"
        "       Then make sure `python3 -V` shows 3.11+ before running\n"
        "       `python3 InstallAll.py ...` again.\n"
        .format(*sys.version_info[:3])
    )
    sys.exit(1)

import shutil
import pathlib
import time
import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(
    prog='InstallAll', 
    description=(
        "Install ecalj and run tests.\n"
        "This script will build ecalj, install binaries and scripts to the specified directory, and run installation tests.\n"
        "Example usage:\n"
        "InstallAll.py --fc ifx --clean --bindir ~/bin  <-- CPU version\n"
        "InstallAll.py --fc nvfortran --clean --bindir ~/bin --gpu --gemmul8 <-- GPU version \n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

parser.add_argument("-np", help='number of mpi cores for install test', default=8, type=int)
parser.add_argument('--clean', help='Clean CMakeCache CMakeFiles before make', action='store_true')
parser.add_argument('--gpu', help='nvfortran for GPU', action='store_true')
parser.add_argument('--bindir', help='ecalj binaries and scripts', type=str, default=str(Path.home() / 'bin'))
parser.add_argument('--fc', help='fortran compiler gfortran/ifort/ifx/nvfortran', type=str, required=True)
parser.add_argument('--notest', help='no test. only compile', action='store_true')
parser.add_argument('--no-bashrc', help='do not append source line to ~/.bashrc', action='store_true')
parser.add_argument('--verbose', help='verbose on for debug', action='store_true')
parser.add_argument('--debug', help='debug', action='store_true')
parser.add_argument('--mp', help='Use mixed precision for test', action='store_true')
parser.add_argument('-np2', help='MPI size for GPU GW executables (default: same as -np)', default=None, type=int)
parser.add_argument('--gemmul8', help='build and install GEMMul8 library (this option is ignored unless --gpu is set)',
                    action='store_true', default=False)
args = parser.parse_args()
args.gemmul8 = args.gpu and args.gemmul8

def run_shell(command, cwd=None, env=None, skip_on_error=False):
    """Executes a shell command and handles failure based on skip_on_error."""
    try:
        # Path objects are automatically converted to strings for subprocess.
        subprocess.run(command, shell=True, check=True, cwd=cwd, env=env)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e.cmd}", file=sys.stderr)
        if not skip_on_error:
            sys.exit(1)
        else:
            print("Skipping command.", file=sys.stderr)

def build_and_install_gemmul8(build_dir: Path, bin_dir: Path):
    repo_url = "https://github.com/RIKEN-RCCS/GEMMul8"
    clone_dir = build_dir / "GEMMul8"
    libfile = clone_dir / "GEMMul8" / "lib" / "libgemmul8.so"
    if not clone_dir.exists():
        run_shell(f"git clone {repo_url} {clone_dir}", skip_on_error=True)
    if not libfile.is_file():
        run_shell("make -j", cwd=clone_dir / "GEMMul8", skip_on_error=True)
    try:
        shutil.copy(libfile, bin_dir)
    except Exception as e:
        print(f"Warning: Failed to copy {libfile} to {bin_dir}: {e}", file=sys.stderr)

ECALJ_BASHRC_MARKER = "# >>> ecalj bash completion (auto-installed by InstallAll.py) >>>"
ECALJ_BASHRC_END    = "# <<< ecalj bash completion <<<"

# Manifest filename dropped into BIN_DIR.  uninstall.py reads this to
# learn which symlinks / real files InstallAll.py placed in BIN_DIR
# and what `~/.bashrc` it appended to.
ECALJ_MANIFEST_NAME = "ecalj_install_manifest.txt"


def write_install_manifest(bin_dir: Path, ecalj_root: Path):
    """Scan BIN_DIR and record every path that came from this ecalj install.

    Includes:
      * symlinks under bin_dir whose target resolves into ecalj_root
        (covers SRC/exec, StructureTool, GetSyml, ecalj_auto helpers,
        and the cmake `deliver` target's libecaljF*.so + main exes).
      * known real files (ecalj_cmdopts.list, libgemmul8.so*).

    Format: header lines starting with '#', then one absolute path per
    line.  Re-running InstallAll.py overwrites the manifest, so stale
    entries from an older install layout do not accumulate.
    """
    manifest = bin_dir / ECALJ_MANIFEST_NAME
    ecalj_root_str = str(ecalj_root)
    entries = []
    real_file_names = {"ecalj_cmdopts.list", ECALJ_MANIFEST_NAME}
    for entry in sorted(bin_dir.iterdir()):
        if entry.is_symlink():
            try:
                target = os.readlink(entry)
            except OSError:
                continue
            target_abs = (bin_dir / target).resolve() if not os.path.isabs(target) else Path(target)
            if str(target_abs).startswith(ecalj_root_str + os.sep) or str(target_abs) == ecalj_root_str:
                entries.append(str(entry))
        elif entry.is_file():
            if entry.name in real_file_names or entry.name.startswith("libgemmul8."):
                entries.append(str(entry))
    timestamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    bashrc = pathlib.Path.home() / ".bashrc"
    bashrc_touched = bashrc.exists() and ECALJ_BASHRC_MARKER in bashrc.read_text()
    lines = [
        "# ecalj install manifest",
        f"# created: {timestamp}",
        f"# ecalj_root: {ecalj_root}",
        f"# bin_dir: {bin_dir}",
        f"# bashrc: {bashrc} (marker {'present' if bashrc_touched else 'absent'})",
        "# uninstall: run `python3 {ecalj_root}/uninstall.py` to remove everything below.".format(ecalj_root=ecalj_root),
        "",
    ]
    lines.extend(entries)
    manifest.write_text("\n".join(lines) + "\n")
    print(f"Wrote install manifest ({len(entries)} entries) -> {manifest}")


def install_bash_completion(bin_dir):
    """Append a guarded source line to ~/.bashrc so tab-completion for
    Legacy2toml.py, lmf, lmfa, lmchk, gwsc, ... is available in new shells.

    Idempotent: if the marker is already present (any earlier install),
    the bashrc is left untouched.
    """
    bashrc = pathlib.Path.home() / ".bashrc"
    snippet = bin_dir / "ecalj_complete.bash"
    if not snippet.exists():
        print(f"Skipping bash completion: {snippet} not found.")
        return
    if bashrc.exists() and ECALJ_BASHRC_MARKER in bashrc.read_text():
        print(f"Bash completion already registered in {bashrc} (marker found).")
        return
    block = (
        f"\n{ECALJ_BASHRC_MARKER}\n"
        "# Tab-complete <sname> for ecalj scripts based on cwd contents.\n"
        "# Remove this block (and the matching end marker) to disable.\n"
        f"[ -f {snippet} ] && source {snippet}\n"
        f"{ECALJ_BASHRC_END}\n"
    )
    with open(bashrc, "a") as f:
        f.write(block)
    print(f"Appended ecalj bash-completion source line to {bashrc}")
    print(f"  -> open a NEW shell, or run:  source {snippet}")


def main():
    if args.gpu:
        import fcntl
        try:
            lockfile = open('/tmp/gpu.lock', 'w')
        except PermissionError:
            print("ERROR: cannot open /tmp/gpu.lock for writing (probably owned by another user).")
            sys.exit(1)
        try:
            fcntl.flock(lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("ERROR: GPU is locked by another job (see /tmp/gpu.lock). Wait or kill the other job.")
            sys.exit(1)

    BUILD_TYPE = "Debug" if args.debug else "Release"
    CWD = Path.cwd()
    BIN_DIR = Path(args.bindir).expanduser().resolve()
    ncore = args.np
    FC = args.fc
    verbose = 'VERBOSE=1 ' if args.verbose else ''

    BIN_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Going to install required binaries and scripts to {BIN_DIR}")
    start_time = time.time()

    # --- Directories ---
    EXEC_DIR = CWD / 'SRC' / 'exec'   # workflow scripts (job_xxx etc.) + pylib
    BUILD_DIR   = CWD / 'SRC' / f'build_{FC}'

    # --- Symlink everything in exec/ into ~/bin ---
    for item in EXEC_DIR.iterdir():
        link = BIN_DIR / item.name
        if link.is_symlink() or (link.exists() and link.is_file()):
            link.unlink()
        elif link.is_dir() and not link.is_symlink():
            shutil.rmtree(link)
        link.symlink_to(item.resolve())

    # --- Extra symlinks into ~/bin from other repo dirs ---
    scripts_to_link = ['StructureTool/viewvesta', 'StructureTool/ctrl2vasp', 'StructureTool/vasp2ctrl', 'GetSyml/getsyml']
    for scr_path_str in scripts_to_link:
        script_path = Path(scr_path_str)
        src_file = CWD / f"{scr_path_str}.py"
        link_path = BIN_DIR / script_path.name
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(src_file)
    print(f"Linked scripts into {BIN_DIR}")

    # ecalj_auto/ helper scripts: keep .py extension in BIN_DIR symlink so callers can
    # `~/bin/slot_run.py` directly (matches how worker.sh / run_gw1500.sh invoke them).
    ecalj_auto_scripts = ['slot_run.py', 'slot_scheduler_daemon.py']
    auto_dir = CWD / 'ecalj_auto'
    for fname in ecalj_auto_scripts:
        src_file = auto_dir / fname
        if not src_file.exists():
            print(f"Warning: ecalj_auto/{fname} not found, skipping symlink", file=sys.stderr)
            continue
        link_path = BIN_DIR / fname
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(src_file)
    print(f"ecalj_auto helper symlinks created in {BIN_DIR}")

    # --- Clean up build directory if requested ---
    if args.clean:
        print("Cleaning previous build files...")
        shutil.rmtree(BUILD_DIR, ignore_errors=True)

    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # --- Configure and Build using CMake ---
    cmake_env = os.environ.copy()
    cmake_env['FC'] = FC

    # Pass BIN_DIR to CMake so the `deliver` target auto-deploys
    # libecaljF*.so + every main exe to BIN_DIR on every build, atomically.
    cmake_options = (
        f"-S {CWD / 'SRC'} -B {BUILD_DIR}"
        f" -DCMAKE_BUILD_TYPE={BUILD_TYPE}"
        f" -DECALJ_BIN_DIR={BIN_DIR}"
    )
    if args.gemmul8:
        build_and_install_gemmul8(BUILD_DIR, BIN_DIR)
    if args.gpu:
        print("Configuring for GPU build...")
        cmake_options += " -DBUILD_MP=ON -DBUILD_GPU=ON -DBUILD_MP_GPU=ON"

    run_shell(f"cmake {cmake_options}", env=cmake_env)

    jobs = min(os.cpu_count(), 8)  # nvfortran ICE with high parallelism
    print(f"Building with {jobs} parallel jobs...")

    run_shell(f"{verbose}cmake --build {BUILD_DIR} -j{jobs}", env=cmake_env)

    # Dump the cmdopt registry into BINDIR so ecalj_complete.bash can
    # offer flag-name completion without exec()ing lmf on every tab.
    # Regenerated on every InstallAll.py run, so a registry edit is
    # picked up the next time the user reinstalls.
    #
    # Wrap in `mpirun -np 1` because some MPI flavours (e.g. HPCX on
    # kt1) refuse to bring up an MPI binary without launcher framing,
    # and abort before it can reach c0_listcmdopt.
    cmdopt_list = BIN_DIR / 'ecalj_cmdopts.list'
    print(f'Dumping cmdopt registry -> {cmdopt_list}')
    # skip_on_error=True: some MPI launchers (HPCX on kt1) flag a non-zero
    # exit even when MPI_Init/MPI_Finalize bracket the dump and stdout is
    # fully written. Tab completion is cosmetic; never let it block the
    # install.
    run_shell(f"mpirun -np 1 {BIN_DIR / 'lmf'} --listcmdopt > {cmdopt_list}",
              skip_on_error=True)

    # Install per-user bash completion (one-shot append to ~/.bashrc).
    if not args.no_bashrc:
        install_bash_completion(BIN_DIR)

    # Record every BIN_DIR entry that belongs to this install so
    # uninstall.py can undo it without re-deriving the layout.
    write_install_manifest(BIN_DIR, CWD)

    if args.notest:
        print('Compilation finished. Skipping tests.')
        return

    # --- Run Install test ---
    print('\n=== Running installation test ===')
    test_dir = CWD / 'Samples' / 'TestInstall'
    end_time_make = time.time()
    start_time_test = time.time()

    test_opts = f"-np {ncore} --all"
    if args.gpu: test_opts += " --gpu"
    if args.mp:  test_opts += " --mp"
    if args.np2: test_opts += f" -np2 {args.np2}"
    run_shell(f"{BIN_DIR / 'testecalj'} {test_opts}", cwd=test_dir)

    end_time = time.time()
    elapsed_time_make = end_time_make - start_time
    elapsed_time_test = end_time - start_time_test
    elapsed_time = end_time - start_time
    print(f"\nElapsed time for make        : {elapsed_time_make:.0f} seconds")
    print(f"Elapsed time for testecalj.py: {elapsed_time_test:.0f} seconds")
    print(f"Total elapsed time           : {elapsed_time:.0f} seconds")

if __name__ == "__main__":
    main()
