# ecalj のアンインストール {#ecalj-uninstall}

`InstallAll.py` は、インストールしたものの一覧（マニフェスト、`ecalj_install_manifest.txt`）を `BIN_DIR`（既定は `~/bin/`）に書きます。
アンインストールは次の 3 つのコマンドで済みます。

```bash
# 1. Remove BIN_DIR symlinks + the bash-completion block in ~/.bashrc.
cd <path-to-ecalj-checkout>
python3 uninstall.py

# 2. Open a new shell so the completion handlers drop out of the session.
#    (Or: `source ~/.bashrc` in this shell.)

# 3. Delete the source checkout itself.
rm -rf <path-to-ecalj-checkout>
```

これだけです。インストールでは `/usr/local/` や `/etc/` などのシステムの場所に何も書かないので、どの段階でも `sudo` は要りません。

## `uninstall.py` がすること {#what-uninstall-py-actually-does}

`BIN_DIR/ecalj_install_manifest.txt` を読み、そこに書かれたパスを全部消します。消すのは次のものです。

* `InstallAll.py` と CMake の `deliver` ターゲットが `BIN_DIR` に置いた **symlink**: `lmf`、`lmfa`、`lmchk`、`gwsc`、`hsfp0`、
  `hx0fp0`、`hvccfp0`、`libecaljF*.so`、`viewvesta`、`ctrl2vasp`、`vasp2ctrl`、`getsyml`、`slot_run.py`、`slot_scheduler_daemon.py`、
  すべての `job_*`、`ecalj_complete.bash` など。
* 一緒に置いた**実体のファイル**: `ecalj_cmdopts.list` と `libgemmul8.*`（後者は `--gemmul8` を付けてインストールしたときだけ）。

次に、`~/.bashrc` に自動で足した bash の補完のブロック（`# >>> ecalj bash completion (auto-installed by InstallAll.py) >>>` と
`# <<< ecalj bash completion <<<` で囲まれた部分）を取り除き、最後にマニフェストのファイルそのものを消します。

よく使うオプション:

```bash
python3 uninstall.py --dry-run        # show what would be removed, change nothing
python3 uninstall.py -y               # skip the y/N confirmation
python3 uninstall.py --bindir ~/opt   # if you installed to a non-default BIN_DIR
python3 uninstall.py --no-bashrc      # leave ~/.bashrc untouched
```

ソースのチェックアウト（`ecalj/` ディレクトリそのもの。`SRC/build_*/` や `Samples/TestInstall/*_work/` も含む）には `uninstall.py` は
**手を付けません**。これを消すのは、上の 3 番目のコマンド（`rm -rf <ecalj-checkout>`）です。

## マニフェストが無いとき {#if-the-manifest-is-missing}

マニフェストを書くようになる前の古いインストールか、マニフェストを手で消した場合です。やり方は 2 つあります。

1. **InstallAll.py を一度走らせ直して**マニフェストを作り直し、それから `uninstall.py` を走らせる。
2. **手で消す**:

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
