# ecalj/ecalj_auto

This is for automatic batch calculations across many POSCAR files.
`ecalj_auto` was used for generating [ecaljdatabase](https://github.com/tkotani/DOSnpSupplement/blob/main/bandpng.md#band-structure--total-dos).

For the GW1500 production run details and slot scheduler, see:
- [ecalj_auto/README.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README.md) — main usage guide
- [ecalj_auto/README_slot_scheduler.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README_slot_scheduler.md) — TOML input + slot scheduler (the combined GW program `hgw` is called `hgw_combined` in that file)

> ⚠️ **Input** — the input of a job is `ctrlg.<sname>.toml`, written by `ctrlgenToml.py`. Of the drivers in `ecalj_auto`,
> `gw1500_rerun.sh` (POSCAR → `vasp2ctrl` → `ctrlgenToml.py --ssig=0.8` → `gwscconv --gpu --prec=fp32`) runs with the present programs.
> `jobsubmit.py` (`auto/creplot.py`: `ctrlgenM1.py`, `GWinput`, `-vnit=...`) and `worker.sh` / `run_gw1500*.sh`
> (`-v[ham.scaledsigma]=0.8`) still carry the retired syntax, on which the present programs stop; write the overrides as
> `--ctrlg:<section.key>=<value>` before using them. Templates: `jobtemplate{,.kugui,.ohtaka,.ucgw}` for SLURM/PBS dispatch. See [TOML migration](../manual/toml_migration.md).
>
> ⚠️ **GW1500 batch: fp32** — the batch runs `gwscconv --gpu --prec=fp32 --conv-tol 0.1`
> (`--gpu --mp --fp32` is the older spelling of the same precision).
> fp32 では、悪条件な誘電行列 (重元素 + 分子アニオン NO3 / N3 /
> ClO など) でも QSGW が NaN / 発散しない。詳細は
> [gwsc § `--fp32`](../manual/gwsc.md#fp32-2026-06) および
> [GPU マニュアル § 混合精度](../manual/ecaljgpu.md#混合精度-mp-と-fp32)。
> 失敗事例の再現は [`Samples/mptf32problem/`](https://github.com/tkotani/ecalj/tree/main/Samples/mptf32problem)。



### python のImportError が発生する場合

python の version や環境によって以下のようなエラーがでる場合がある(ISSP python3.6で発生)
```text OUTPUT/testSGA/start@xxxxxxx/job0.out
Traceback (most recent call last):
  File "/home/k0413/k041300/ecalj/ecalj_auto/OUTPUT/testSGA/start@20250410-095838/job_mp.py", line 3, in <module>
    import pandas as pd
  File "/home/k0413/k041300/.local/lib/python3.6/site-packages/pandas/__init__.py", line 17, in <module>
    "Unable to import required dependencies:\n" + "\n".join(missing_dependencies)
ImportError: Unable to import required dependencies:
pytz: No module named 'pytz'
```
この場合pythonのversion upを試す。

#### `mise` を使用する場合.　(pyenvでも同様です)
以下を `~/.bashrc` に記載 (zsh の場合は 2行目のbash を zsh とする)
```bash
type mise > /dev/null 2>&1 || curl https://mise.run | sh
eval "$(~/.local/bin/mise activate bash)"
```
>  [!NOTE]
> mise は パッケージ管理ソフトの一種であり、詳細は [mise](https://mise.jdx.dev/) を参考にして下さい。

その後 ~/.bashrc  の再読み込みを行うと mise がインストールされ使用できるようになる。
```bash
source ~/.bashrc
```

ecalj_auto があるディレクトリで以下を実行し, python をinstall する。
```bash 
mise use python@latest
```
> [!IMPORTANT]
> ここでinstallされたpythonは, ecalj_auto があるディレクトリより下位のディレクトリでのみ有効となることに注意。

ここでinstallしたpythonが使用するlibraryを以下で導入する。
```bash
pip3 install pandas seekpath spglib --user
```

## サーバーのコマンドに合わせてスクリプトの変更

### mpirun 関係
MPI の実行コマンドをサーバの使用に合わせて変更する．デフォルトは`mpirun` であり, 以下のスクリプトに記載されている．
ecalj/ecalj_auto/auto/creplot.py
```python
    def run_lmf(self, fout,foute):
        command1 = ['mpirun', '-np', self.ncore, self.epath/'lmf', self.num] + self.option_lmf
        command2 = ['tail', '-f', f'save.{self.num}']
        run_popen(command1, command2, fout, foute, 'a')
        return check_save(f'save.{self.num}')
```
`command1 = ['mpirun', '-np', self.ncore, self.epath/'lmf', self.num] + self.option_lmf` を修正する．
#### SLURM の場合: 例 ISSP system B

On the slurm case: e.g., ISSP system B, Othtaka
```python
        command1 = ['srun', '-n', self.ncore, self.epath/'lmf', self.num] + self.option_lmf
```

#### OpenMPI の場合: 例 ISSP system C
```python
        command1 = ['mpiexec', '--bind-to none', '-np', self.ncore, self.epath/'lmf', self.num] + self.option_lmf
```

### ジョブ投入関係
ジョブ投入コマンドはサーバに合わせて変更する．現在のデフォルトは
**バックグラウンドの `bash`** で, 以下のスクリプトに記載されている．
`~/ecalj/ecalj_auto/auto/Job.py`:

```python
            os.system(f'bash {jobx} &')
```

(2025 後期に `qsub` から `bash` 起動へ migrate された。SLURM/PBS
で投入したい場合は対象行を以下に書き換える。)

#### qsub の場合 (PBS):
```python
            os.system(f'qsub {jobx}')
```

#### SLURM の場合 (例 ISSP system B):
```python
            os.system(f'sbatch {jobx}')
```

## Query

MPからのPROCARの取得方法

依存関係: 以下のpython ライブラリが必要
- pymatgen
- mp-api
インストールはpipから
```bash
pip3 install pymatgen mp-api --user
```
Materials Project の API key が必要。個人のものなので、ecalj の最上位に `MaterialProject.key`（1 行、git には入らない）として置く。
公開用の見本 `MaterialProject.key.example` を写して、置き換え用の行を自分のキーにする。環境変数 `MP_API_KEY` があればそちらが優先。
[`ecalj_auto/config.ini`](../../ecalj_auto/config.ini) はリポジトリに入っているので、キーを書かない（2026-10-01 まではここに書く形で、キーが公開リポジトリに入っていた）。

```bash
cd ~/ecalj; cp MaterialProject.key.example MaterialProject.key; chmod 600 MaterialProject.key   # そして最後の行を自分のキーに
```
