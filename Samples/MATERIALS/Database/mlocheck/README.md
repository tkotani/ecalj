# MLO の自動の模型の試験（2026-10-01）の道具

`Samples/MATERIALS` の全物質で LDA のあと `job_mlo` の既定の模型を作り、対称線の上で DFT のバンドと比べた（結果は `../MLOcheck_20261001.tsv`、
`../MLOcheck_20261001_variants.tsv`、経緯は `MD/research_log.md` の 2026-10-01 朝 07:16）。ここにあるのは、そのときに回したままのスクリプト。

| ファイル | 何をするか |
|---|---|
| `prep.py` | 作業場所（このファイルを置いた場所）に物質ごとのディレクトリを作り、各 ctrlg を写して `[mlo]` を書き直す（`mlo_lm` の規則は冒頭の注釈）。`materials.json` を書く |
| `run_one.sh <物質> <np>` | `lmfa` → `lmf` → `getsyml --nobzview` → `job_band` → `job_mlo` → `mlo_bandplot.py`。`../status.log` に時刻を書く |
| `run_soc.sh <物質> <np>` | so = 1 の物質。最後が `job_mlo_soc`（その出力は `band_MLO_spin1.dat`。スクリプトが探す `.soc.dat` は誤りで、FAIL と書くが結果はできている） |
| `worker.sh <np>` | `queue.txt` から物質を 1 つずつ取って `run_one.sh` を回す（`flock`） |
| `mkvar.py <元> <名前> <種類>` | 終わった物質を写し、`mlo_lm2`（EH2 の s,p）か `mlo_lm3`（半内殻の局所軌道）を足して `job_mlo` だけ回す。種類は `lm2sp`・`lm3d`・`lm3semi`・`both` |
| `gallery.py <出力>` と `page_template.html` | `mlo_bandcheck.py --json` の結果と図から一覧のページ（`index.html` と `img/`）を作る |

評価は `SRC/exec/mlo_bandcheck.py <dir> ... --json out.json --tsv out.tsv`。
