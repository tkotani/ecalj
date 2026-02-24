import sys
import glob
import shutil
import subprocess
from pathlib import Path

def run_shell(command: str, cwd=None, env=None, skip_on_error: bool = False):
    """Executes a shell command and handles failure based on skip_on_error."""
    try:
        subprocess.run(command, shell=True, check=True, cwd=cwd, env=env)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e.cmd}", file=sys.stderr)
        if not skip_on_error:
            sys.exit(1)
        else:
            print("Skipping command.", file=sys.stderr)


def merge_files(pattern: str, output_file: str | Path, remove_sources: bool = True):
    """Merges files matching a pattern into a single output file."""
    output = Path(output_file).resolve()
    files = sorted(Path(p).resolve() for p in glob.glob(pattern))
    files = [f for f in files if f != output]
    if not files:
        return
    with output.open("wb") as wfd:
        for f in files:
            with f.open("rb") as fd:
                shutil.copyfileobj(fd, wfd)
            if remove_sources:
                f.unlink()

def ensure_dir(dirname: str | Path):
    """Ensures a directory exists, creating it if necessary."""
    Path(dirname).mkdir(exist_ok=True)

def remove_files(*patterns: str):
    """Removes files matching given glob patterns."""
    for pattern in patterns:
        for f in sorted(Path().glob(pattern)):
            print(f"Removing {f}")
            f.unlink()

def _resolve_files(*sources: str | Path) -> set[Path]:
    """Helper to expand file paths and glob patterns into a set of file paths."""
    files = set()
    for src in sources:
        files.update(Path().glob(str(src)))
    return files

def copy_files(*sources: str, dest_dir: str | Path):
    """Copies files to a destination directory.

    Sources can be explicit file paths or glob patterns.
    """
    destination = Path(dest_dir)
    destination.mkdir(exist_ok=True)
    
    resolved_files = _resolve_files(*sources)
    for f_path in resolved_files:
        if f_path.is_file():
            shutil.copy(f_path, destination)

def move_files(*sources: str, dest_dir: str | Path):
    """Moves files to a destination directory.

    Sources can be explicit file paths or glob patterns.
    """
    destination = Path(dest_dir)
    destination.mkdir(exist_ok=True)

    resolved_files = _resolve_files(*sources)
    for f_path in resolved_files:
        if f_path.is_file():
            f_path.replace(destination / f_path.name)
