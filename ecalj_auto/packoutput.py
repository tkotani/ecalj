#!/usr/bin/env python3
import os
import sys
import shutil
import fnmatch
from pathlib import Path

def match_patterns(filename, patterns, exclude_patterns):
    return any(fnmatch.fnmatch(filename, pattern) for pattern in patterns) and not any(fnmatch.fnmatch(filename, pattern) for pattern in exclude_patterns)

patterns = ["QPU*", "QPD*", "bandplot*", "bnd*", "ctrl.*", "ctrls.*", "GWinput", "job*", "log*", "save.*", "SiteInfo.*","dos*"]
exclude_patterns = ["*.tmp*", "*temp*", "ctrls.POSCAR*"]
dpatterns = ["mp-*", "QSGW.*", "LDA"]
exclude_dpatterns = []  

readdir = sys.argv[1]
d = Path(readdir)
full_path = d.resolve()
newdir = full_path.parent / (full_path.stem + "_pack")
newdir.mkdir(parents=True, exist_ok=True)

for root, dirs, files in os.walk(readdir):
    rel_path = Path(root).relative_to(readdir)
    new_root = newdir / rel_path
    for dir in dirs:
        if match_patterns(dir, dpatterns, exclude_dpatterns):
            dest_dir = new_root / dir
            print('--- Generate ', dest_dir)
            os.makedirs(dest_dir, exist_ok=True)
    for file in files:
        src_path = Path(root) / file
        if match_patterns(file, patterns, exclude_patterns):
            dest_path = new_root / file
            print(src_path, dest_path)
            shutil.copy2(src_path, dest_path)
