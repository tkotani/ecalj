# ecalj Uninstall

`InstallAll.py` writes an install manifest into `BIN_DIR`
(`ecalj_install_manifest.txt`, default `~/bin/`).  Uninstall is then
three commands:

```bash
# 1. Remove BIN_DIR symlinks + the bash-completion block in ~/.bashrc.
cd <path-to-ecalj-checkout>
python3 uninstall.py

# 2. Open a new shell so the completion handlers drop out of the session.
#    (Or: `source ~/.bashrc` in this shell.)

# 3. Delete the source checkout itself.
rm -rf <path-to-ecalj-checkout>
```

That's it.  Nothing is written to `/usr/local/`, `/etc/`, or any system
path during install, so no `sudo` is needed at any step.

## What `uninstall.py` actually does

It reads `BIN_DIR/ecalj_install_manifest.txt` and removes every path
listed there, which covers:

* **Symlinks** placed in `BIN_DIR` by `InstallAll.py` and by the
  CMake `deliver` target — `lmf`, `lmfa`, `lmchk`, `gwsc`, `hsfp0`,
  `hx0fp0`, `hvccfp0`, `libecaljF*.so`, `viewvesta`, `ctrl2vasp`,
  `vasp2ctrl`, `getsyml`, `slot_run.py`, `slot_scheduler_daemon.py`,
  every `job_*` workflow, `ecalj_complete.bash`, ...
* **Real files** dropped alongside them — `ecalj_cmdopts.list` and
  `libgemmul8.*` (the latter only if you installed with `--gemmul8`).

Then it strips the auto-installed bash-completion block from
`~/.bashrc` (the one bracketed by
`# >>> ecalj bash completion (auto-installed by InstallAll.py) >>>` and
`# <<< ecalj bash completion <<<`), and finally deletes the manifest
file itself.

Common flags:

```bash
python3 uninstall.py --dry-run        # show what would be removed, change nothing
python3 uninstall.py -y               # skip the y/N confirmation
python3 uninstall.py --bindir ~/opt   # if you installed to a non-default BIN_DIR
python3 uninstall.py --no-bashrc      # leave ~/.bashrc untouched
```

The source checkout (the `ecalj/` directory itself, including
`SRC/build_*/` and `Samples/TestInstall/*_work/`) is **not** touched by
`uninstall.py` — step 3 above (`rm -rf <ecalj-checkout>`) is the one
command that removes it.

## If the manifest is missing

You may have an older install that pre-dates the manifest (or you
removed the manifest by hand).  Two options:

1. **Re-run InstallAll.py once** to regenerate the manifest, then run
   `uninstall.py`.
2. **Remove by hand**:

   ```bash
   ECALJ_DIR=~/ecalj     # path to your ecalj checkout
   BIN_DIR=~/bin         # whatever you passed as --bindir

   # All BIN_DIR entries from ecalj are symlinks back into the clone.
   find "$BIN_DIR" -maxdepth 1 -type l -lname "${ECALJ_DIR%/}/*" -delete
   rm -f "$BIN_DIR/ecalj_cmdopts.list" "$BIN_DIR"/libgemmul8.*

   # Strip the bash-completion block.
   sed -i '/# >>> ecalj bash completion (auto-installed by InstallAll.py) >>>/,/# <<< ecalj bash completion <<</d' ~/.bashrc

   # Remove the checkout.
   rm -rf "$ECALJ_DIR"
   ```
