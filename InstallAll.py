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
    SCRIPTS_DIR = CWD / 'SRC' / 'scripts'   # Python/shell scripts source
    EXEC_DIR    = CWD / 'SRC' / 'exec'      # assembled runtime (scripts + binary symlinks)
    BUILD_DIR   = CWD / 'SRC' / f'build_{FC}'

    # --- Assemble exec/: symlink every file/dir in scripts/ into exec/ ---
    EXEC_DIR.mkdir(parents=True, exist_ok=True)
    for item in SCRIPTS_DIR.iterdir():
        link = EXEC_DIR / item.name
        if link.is_symlink() or (link.exists() and link.is_file()):
            link.unlink()
        elif link.is_dir() and not link.is_symlink():
            shutil.rmtree(link)
        link.symlink_to(item.resolve())

    # --- Extra symlinks into exec/ and ~/bin from other repo dirs ---
    scripts_to_link = ['StructureTool/viewvesta', 'StructureTool/ctrl2vasp', 'StructureTool/vasp2ctrl', 'GetSyml/getsyml']
    for scr_path_str in scripts_to_link:
        script_path = Path(scr_path_str)
        src_file = CWD / f"{scr_path_str}.py"
        for dest_dir in [BIN_DIR, EXEC_DIR]:
            link_path = dest_dir / script_path.name
            if link_path.exists() or link_path.is_symlink():
                link_path.unlink()
            link_path.symlink_to(src_file)
    print(f"Assembled exec/ and created links in {BIN_DIR}")

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

    # --- Copy scripts from SCRIPTS_DIR to BIN_DIR ---
    # libecaljF*.so + every main binary are deployed atomically by cmake deliver.
    # Here we copy the Python/shell scripts from scripts/.
    print(f'Copying helper scripts to {BIN_DIR}')
    for path_item in SCRIPTS_DIR.iterdir():
        if path_item.name == 'pylib':
            continue
        if path_item.is_file() and path_item.suffix != '.so' and os.access(path_item, os.X_OK):
            try:
                shutil.copy2(path_item, BIN_DIR)
            except (OSError, PermissionError) as e:
                print(f"Warning: Skipping {path_item.name}: {e}", file=sys.stderr)

    # Copy clusters.toml
    for data_file in ('clusters.toml', 'ecalj_complete.bash'):
        src = SCRIPTS_DIR / data_file
        if src.exists():
            shutil.copy2(src, BIN_DIR)
            print(f"Copied {data_file} to {BIN_DIR}")

    # Copy pylib/ package
    pylib_src = SCRIPTS_DIR / 'pylib'
    if pylib_src.is_dir():
        shutil.copytree(pylib_src, BIN_DIR / 'pylib', dirs_exist_ok=True)
        print(f'Copied pylib/ to {BIN_DIR}')

    # Install per-user bash completion (one-shot append to ~/.bashrc).
    if not args.no_bashrc:
        install_bash_completion(BIN_DIR)

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
