#!/usr/bin/env python3
import os
import sys
import shutil
import time
import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(prog='InstallAll', description='''
Install ecalj and run tests.
''')
parser.add_argument("-np", help='number of mpi cores for install test', default=8, type=int)
parser.add_argument('--clean', help='Clean CMakeCache CMakeFiles before make', action='store_true')
parser.add_argument('--gpu', help='nvfortran for GPU', action='store_true')
parser.add_argument('--bindir', help='ecalj binaries and scripts', type=str, default=str(Path.home() / 'bin'))
parser.add_argument('--fc', help='fortran compiler gfortran/ifort/ifx/nvfortran', type=str, required=True)
parser.add_argument('--notest', help='no test. only compile', action='store_true')
parser.add_argument('--verbose', help='verbose on for debug', action='store_true')
parser.add_argument('--debug', help='debug', action='store_true')
args = parser.parse_args()

def run_command(command, cwd=None, env=None, skip_on_error=False):
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
        run_command(f"git clone {repo_url} {clone_dir}", skip_on_error=True)
    if not libfile.is_file():
        run_command("make -j", cwd=clone_dir / "GEMMul8", skip_on_error=True)
    try:
        shutil.copy(libfile, bin_dir)
    except Exception as e:
        print(f"Warning: Failed to copy {libfile} to {bin_dir}: {e}", file=sys.stderr)

def main():
    BUILD_TYPE = "Debug" if args.debug else "Release"
    CWD = Path.cwd()
    BIN_DIR = Path(args.bindir).expanduser().resolve()
    ncore = args.np
    FC = args.fc
    verbose = 'VERBOSE=1 ' if args.verbose else ''

    BIN_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Going to install required binaries and scripts to {BIN_DIR}")
    start_time = time.time()

    # --- Make links ---
    EXEC_DIR = CWD / 'SRC' / 'exec'
    BUILD_DIR = EXEC_DIR / 'build'

    scripts_to_link = ['StructureTool/viewvesta', 'StructureTool/ctrl2vasp', 'StructureTool/vasp2ctrl', 'GetSyml/getsyml']
    for scr_path_str in scripts_to_link:
        script_path = Path(scr_path_str)
        src_file = CWD / f"{scr_path_str}.py"

        for dest_dir in [BIN_DIR, EXEC_DIR]:
            link_path = dest_dir / script_path.name
            if link_path.exists():
                link_path.unlink()
            link_path.symlink_to(src_file)
    print(f"Symbolic links created in {BIN_DIR} and {EXEC_DIR}")

    # --- Clean up build directories if requested ---
    if args.clean:
        print("Cleaning previous build files...")
        makefile = EXEC_DIR / 'Makefile'
        if makefile.exists():
            run_command("make clean", cwd=EXEC_DIR)

        (EXEC_DIR / 'CMakeCache.txt').unlink(missing_ok=True)
        shutil.rmtree(EXEC_DIR / 'CMakeFiles', ignore_errors=True)
        shutil.rmtree(BUILD_DIR, ignore_errors=True)

    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # --- Configure and Build using CMake ---
    cmake_env = os.environ.copy()
    cmake_env['FC'] = FC

    cmake_options = f"-S {EXEC_DIR} -B {BUILD_DIR} -DCMAKE_BUILD_TYPE={BUILD_TYPE}"
    if args.gpu:
        print("Configuring for GPU build...")
        build_and_install_gemmul8(BUILD_DIR, BIN_DIR)
        cmake_options += " -DBUILD_MP=ON -DBUILD_GPU=ON -DBUILD_MP_GPU=ON"

    run_command(f"cmake {cmake_options}", env=cmake_env)

    jobs = min(os.cpu_count(), 32)
    print(f"Building with {jobs} parallel jobs...")
    run_command(f"{verbose}cmake --build {BUILD_DIR} -j{jobs}")

    # --- Copy executables to BIN_DIR ---
    print(f'Copying executables to {BIN_DIR}')
    for d in (EXEC_DIR, BUILD_DIR):
        for path_item in d.iterdir():
            if path_item.is_file() and os.access(path_item, os.X_OK):
                try:
                    shutil.copy(path_item, BIN_DIR)
                except (OSError, PermissionError) as e:
                    print(f"Warning: Skipping {path_item.name}: {e}", file=sys.stderr)

    # Copy clusters.toml from EXEC_DIR to BIN_DIR
    clusters_toml_src = EXEC_DIR / 'clusters.toml'
    if clusters_toml_src.exists():
        try:
            shutil.copy(clusters_toml_src, BIN_DIR)
            print(f"Copied {clusters_toml_src} to {BIN_DIR}")
        except (OSError, PermissionError) as e:
            print(f"Warning: Failed to copy {clusters_toml_src} to {BIN_DIR}: {e}", file=sys.stderr)
    else:
        print(f"Info: {clusters_toml_src} not found, skipping copy.")

    if args.notest:
        print('Compilation finished. Skipping tests.')
        return

    # --- Run Install test ---
    print('\n=== Running installation test ===')
    test_dir = CWD / 'Samples' / 'TestInstall'
    end_time_make = time.time()
    start_time_test = time.time()

    run_command(f"{BIN_DIR / 'testecalj'} -np {ncore} --all", cwd=test_dir)

    end_time = time.time()
    elapsed_time_make = end_time_make - start_time
    elapsed_time_test = end_time - start_time_test
    elapsed_time = end_time - start_time
    print(f"\nElapsed time for make        : {elapsed_time_make:.0f} seconds")
    print(f"Elapsed time for testecalj.py: {elapsed_time_test:.0f} seconds")
    print(f"Total elapsed time           : {elapsed_time:.0f} seconds")

if __name__ == "__main__":
    main()
