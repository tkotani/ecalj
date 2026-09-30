# ecalj 過去ログ — 片付けたものに含まれていたノウハウと、役に立つかもしれない情報

リポジトリから外したもの（各計算機の `ecalj/trash/` に移した）の中にあったノウハウや経緯を、話題ごとにまとめる。
入口は [CLAUDE.md](../CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）。
利用者向けの説明は ecaljdoc が軸なので、ここには開発の側の細かいことだけを書く。

外したファイルの中身は git の履歴から取り出せる（表 1 の「外す前のコミット」を使う）:

```bash
git show <コミット>:<パス>                 # 1 ファイルを見る
git checkout <コミット> -- <パス>          # 作業ツリーに戻す（戻したら追跡し直すかを決める）
```

`trash/` は `.gitignore` に入っていて、空にする（本当に消す）のはメンテナが決める。

---

## 1. ビルドとインストール

### 1.1 `InstallAll.py` の打ち込み（最上位の `g`・`gg`・`i`・`ii`・`n`、2025-10〜2026-02）

1 行だけのシェルの控えだった。いまの使い方は ecaljdoc の ForDevelopers §4（2026-05 の流れは MD/README.md）。

| ファイル | 中身 |
| --- | --- |
| `g` | `./InstallAll.py --fc gfortran` |
| `gg` | `./InstallAll.py --fc gfortran --notest` |
| `i` | `./InstallAll.py --fc ifx` |
| `ii` | `./InstallAll.py --fc ifort` |
| `n` | `./InstallAll.py --fc nvfortran --gpu` |

### 1.2 nvfortran が `gaugm.f90` で内部エラーになった回避（`build_nvfortran.sh`、2026-03-25）

- 当時の nvfortran は `gaugm.f90` のコンパイルで内部エラー（ICE）を出した。回避は、4 つの版（無印・`-D__MP`・`-D__GPU`・両方）ごとに
  `mpifort -O1 -fpic -Mbackslash -cpp`（`-acc`・`-Mvect` なし）で `gaugm.f90` だけを手でコンパイルし、`gaugm.f90` のタイムスタンプを
  古くしてから `gmake` を続け、最後にタイムスタンプを戻す、というもの
- 2026-09 の nvfortran 26.1 では、この回避なしで `InstallAll.py` が通る（kt1・kr7 の試験）。`InstallAll.py` にも CMake にも回避は入っていない
- nvfortran の別の症状（ecaljclaude.md「ビルド環境」）: 毎回違うファイルで signal 11 になることがある。同じコマンドをやり直せば通る

### 1.3 ISSP のスーパーコンピュータでのインストール（`jobinstall_kugui.sh`・`jobinstall_ohtaka.sh`）

インストール（コンパイルと試験）をジョブとして投げるもの。フロントエンドでは MPI を走らせられないため。詳しい手順は ecaljdoc の `UsageISSP.md`。
ベンチマークの投入スクリプトは `Samples/BenchmarkTest/ISSP/`。

- kugui（システム C）: `#PBS -q i2cpu`、`select=1:ncpus=16:mpiprocs=8:ompthreads=1`。
  `ulimit -s unlimited; module purge; module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl`、`./InstallAll.py --fc nvfortran --gpu`
  （GPU ノードなら `-q i1accs`、`ncpus=64:mpiprocs=64`）
- ohtaka（システム B）: `#SBATCH -p i8cpu -N 1 -n 8 -c 1 --ntasks-per-node=8`、`./InstallAll.py --fc ifort`

### 1.4 CMake 以前のビルドと古いソースの控え（`SRC/exec/BK`、〜2026-09）

- `Make.inc.{gfortran,gfortran_mpik,ifort,ifort_mpik}`、`Makefile`・`Makefile.issp`・`Makefile.mac`: Makefile でビルドしていた頃の設定
- `*.F`（`asars.F`、`hambls.F`、`m_rdctrl.F`、`potpus.F` など 26 本）: 固定形式の頃のソースの控え。今のソースは `SRC/subroutines/*.f90`
- `gwsc.py`・`gwsc_org`・`gwsc_sym`・`gwscempty`、`run_arg2`・`run_argloc`: `gwsc` と `run_arg` の古い版

### 1.5 SRC の中のビルドの生成物（2026-10-01 に外した）

`SRC/exec/build`・`SRC/exec/build_gf14`（CMake の作業ツリー）、`SRC/exec_gfortran`・`SRC/exec_gfortran-14`・`SRC/exec_ifx`（`libecaljF.so` の写し）の
1853 ファイル・78 MB が、2026-09-24 のコミット `a0c7a7300` で誤って追跡されていた（公開のリポジトリへの push の前に見つけた）。
`.gitignore` に `/SRC/exec/build*/`・`/SRC/exec_*/` を足した。ビルドの置き場は `SRC/build_<コンパイラ>`（既に無視されている）。

## 2. GPU と計算機

### 2.1 GPU のロック（`gpu_wait_run.sh`、2026-04-29）

- `/tmp/gpu.lock` を `flock` で取ってから引数のコマンドを実行する、1 本の GPU を順番に使うための包み
- いまは `SRC/exec/pylib/gpu_lock.py`（GPU ごとに `/tmp/ecalj_res/gpu<N>.lock`）を `pylib/run_cmd.py` が取る。`CUDA_VISIBLE_DEVICES` の GPU だけ、
  `ECALJ_GPU_LOCK=0` で切れる

## 3. リファクタの記録

### 3.1 `m_sxcf_sc` の状態の整理の事前調査（`.refactor_notes/`、2026-04-30）

`m_itq` のリファクタ（`5cebb798`）の後、streaming（Phase 1-C）をやり直すために、`m_sxcf_sc.f90` の module 変数を
どこまで安全に動かせるかを調べたもの。

- OpenACC の directive がかかる変数（`zsecall`、`wvi_upper`・`wvr_upper`、`idx_i`・`idx_j`）は動かすと危ない。
  directive の無い、相関の側だけで使う変数（`omega`、`keepwv`、`nt0p`・`nt0m`、`ifrcw`・`ifrcwi`、ストップウォッチ）は安全。
  交換と相関の両方で使う変数（`ekc`、`eq`、`ntqxx`、`wkkr`）は両方を同時に書き換える必要がある
- ループの添字（`kx`、`irot`、`ip`、`isp`）は、BLOCK の host scope を確かめたうえで、各サブルーチンのローカルにして安全
- 回帰の道具 `snapshot_baseline.sh`: 各 `*_gwsc_work` の `SECU`・`SEC2U`・`SEXU`・`SEX2U` の md5 と、`lcombined` の診断行
  （`ntq`、`nbandmx`、`ncount` など）を段階ごとに記録して比べる（`snapshots/pre_step2`〜`stepWB1`）
- 調査では derived type にまとめる段階（Step 2.1〜2.5）を提案したが、最終的には type を使わず module の singleton にした
  （WB.3、ecaljclaude.md「実例: WB.3 リファクタ」）

### 3.2 `hrcxq` と `hsfp0_sc` の統合、W のメモリ中継（`PHASE1B_REFACTOR.md`、2026-04-30）

`hrcxq`（W）と `hsfp0_sc --job=2`（Σ_c）を 1 つの実行ファイル（当時 `hgw_combined`、いまの `hgw`）にまとめ、間の `__WVR.<iq>`・`__WVI.<iq>`
（1 物質あたり数十 GB）をメモリの中で受け渡すようにした。構成と段階（WB.1〜WB.3e）は ecaljclaude.md「hgw_combined」。
当時の測定: CPU の試験（si・nio・fe・gas）で QPU がビット一致、GW1500（kt1、`--gpu --mp`）で 1 物質 12〜30% 速く、ディスク 46〜87 GB 減。
リファクタの前の状態はタグ `pre-hgw_combined-baseline`。設計の判断は次のとおり（ほかの統合にも使える）:

- **`Genallcf_v3` を 1 回で済ませる**: `hrcxq` は `incwfx = 0`（ForX0 の core）、`hsfp0_sc` は `incwfx = -1`（ForSxc の core）で呼んでいた。
  両方 `-1` に揃えても si・nio・fe・gas で QPU がビット一致だった（試験の GWinput では ForX0 と ForSxc の列が同じだったため）
- **W の受け渡しは module の buffer とフラグ**: optional 引数（`output_target`、`wv_buf`）を `W0w0i` や `sxcf_scz_correlation` まで
  伝えると呼び出しの形が汚れるので、`m_llw` に buffer と `wv_in_memory` を置き、使う側は `use m_llw` で読む
- **1 ランクだけが計算した W の共有は `MPI_Allreduce(SUM)`**: 各 iq の W は 1 つのランクだけが計算し、ほかはその枠が 0。和を取れば全員が持てる。
  `MPI_Allgatherv` より簡単。ノード内なら `MPI_Win_allocate_shared` でさらに減らせる（未着手、TODOandQuestion.md）
- **Fortran の `.and.` は短絡しない**: `if (present(x) .and. x)` は x が無いときも x を読んで落ちる（このリファクタで一度 SIGSEGV）。
  `if (present(x)) then; if (x) ...` と入れ子にするか、`y = .false.; if (present(x)) y = x`
- **module を 2 回使うときの分類**: (A) 引数で大きさが変わるもの（`Mptauof_zmel`、`Rdpp`、`setitq_hsfp0sc`）は解放して確保し直す、
  (B) ファイルを読むだけのもの（`readhamindex0`、`Genallcf_v3`）は「済み」のガードで 2 回目を飛ばす、
  (C) もともと対応済み（`Read_BZDATA`、`set_m2e_prod_basis`）。全部を片付けて読み直すのは無駄
- 精度に応じた MPI の型: `#ifdef __MP` で `MPI_COMPLEX`、そうでなければ `MPI_DOUBLE_COMPLEX`（`m_llw`）

## 4. サンプルと試験

### 4.1 MLOsamples の古い試行（2026-10-01 に外した）

`Al2O3_Cr/test1`〜`test7`、`RuO2/temp`、`Si666gwsc/temp`、`SrTiO3/temp`、`Al2O3_Cr/ctrl.tmp`、`MLOsamples/GWinput.tmp`、各サンプルの
`bandplot.isp*.glt.bk`。旧形式の入力（`GWinput`）で MLO のパラメータを試したときの跡で、試験（`test.py`）からも README からも参照されていなかった。

### 4.2 `Samples/TestInstall/TESTunused`（cdte、er、gaslc、srtio3、tio2、2025-10）

`testecalj` の対象から外した試験の入力。どの試験の道具からも参照されていなかった。

### 4.3 ecaljdoc の古い下書き（`ecaljdoc_drafts/`、2026-06）

`gw_diagnostics.md`・`gw_lmfh_qpe.md`・`tetrakbt.md`。どれも冒頭に「古い草稿。いまの説明は ecaljdoc の manual/kBT・manual/gwsc」とあり、置き換え済み。
中身の独自の部分は §5 に移した。

## 5. 一発 GW（`gw_lmfh`）と任意の k 線の QP エネルギー（`FiniteT_and_QPE_HOWTO.md`、2026-06-13）

2026-06 の利用者向けの手引き。有限温度の χ₀ のキーはその後変わった（いまは `t_tetrakbt` が必須で 0 / T > 0 / −T、ecaljdoc manual/kBT）。
`gw_lmfh` の使い方、`--dumpW`・`--WVR2ptRaxis`・`HistBin_*` は ecaljdoc（gwsc.md、cmdopts.md、kBT.md）にある。ecaljdoc に無いものをここに残す:

- **メッシュの外の q（`QforGW`）では `EMAXforGW` が必須**: 固有値の解法が nband より少ない状態しか返さず、残りを 1d20 で埋める。
  `EMAXforGW` が無いと、この詰め物の状態まで Σ の対象になり `hsfp0` が `sxcf 222: |w-e| out of range` で止まる
- **窓（`EMINforGW`・`EMAXforGW`）を変えたら交換からやり直す**: `SEXU`・`XCU`（`hsfp0 --job=3`・`--job=11`）と `SECU`（`--job=12`）の
  状態の数が合わないと、`hqpe` が `SEXU` の EOF で止まる
- **`--mp` だけ（TF32）で W を作らない**: 条件の悪い誘電行列が壊れる（2026-06-04 の mp-8196）。`--fp32`（本当の FP32 の積）を併せる。
  単精度で書いた `WVR`・`WVI` は、倍精度の `hsfp0` が読むときに昇格させる（`hsfp0: single-precision WVR/WVI detected`）
- **CPU の相関のメモリ**: QPE のモードでは各ランクが W(ω) を全部持つ（約 ngb²·nω·16 B。LiTi₂O₄ 6³ で 1 ランク 4.6 GB）。
  256 GB の機械で 60 ランクはメモリ不足、20 ランクなら入った。GPU 版はデバイスに置く
- **`gw_lmfh` は対角の Σ_nn だけ**: 帯が交わるところの折れは、縮退した部分空間の非対角の Σ から来る。
  2×2 の関係 (E₊ − E₋)² = (d₁ − d₂)² + 4|Σ₁₂|² に、LDA の分裂・対角 QPE の分裂・QSGW のバンドの分裂（どれも基準の要らない差）を入れて |Σ₁₂| を見積もれる
- **小さい q の誘電関数をただで見る**: `hgw`・`hx0fp0` の標準出力に、ずらした Γ の q について `epsWVR: iq iw omg eps(wFC) eps(woLFC)` の行がある。
  grep すればプラズモン極（Re ε = 0）が分かる
- **実軸のメッシュ**: 準位の ω がプラズモンの鋭い構造に落ちる状態は、ビンの幅に敏感。`HistBin_ratio` を 1.03 → 1.015 と細かくして確かめる。
  ほとんど離散的な極は、細かくするだけでは収束しない。先に幅を与える（`t_tetrakbt > 0`、またはビンより広い Gaussian の `t_tetrakbt < 0`）
- W のファイルの名前: `__WV.d`、`__WVR.<iq>`、`__WVI.<iq>`（接頭辞の無い `WV.d` は古い名前）

## 6. 2026-06〜09 の更新の要点（`HIGHLIGHTS_2026-06_09.md`、2026-09-30 時点の要約）

2026-06-03 以降に main に入った 140 コミットの要点。詳細は `Changes.txt`（日付順）と ecaljdoc。

**6.1 入力は `ctrlg.<sname>.toml` 一つになった（破壊的変更）**

- lmf も GW も `ctrlg.<sname>.toml` だけを読む。`PB.<sname>.toml`、`esm_input.dat`、`GWinput.toml` は廃止。前二者が残っていると止まり、
  変換（`ctrlg_absorb.py <sname>`、原本は `.bk`）を案内する。旧 `ctrl.<sname>` + `GWinput` からは `Legacy2toml.py <sname>`
- 節の順は `[gw]` `[mlo]` `[blocks]` `[product_basis]`（末尾に nlx / valence / core の表）。`Worb` → `[mlo] mlo_lm`、`QforEPS`・`QforGW` → `[gw]`。
  空行は `# === X ===` の大見出しの前だけ、節の中は `# ----` の罫線（`ctrlgenToml.py`・`gwinit`・`Legacy2toml.py` が同じ形を出す）
- 新規は `ctrls.<sname>`（構造だけ）→ `ctrlgenToml.py <sname>` → `lmfa`・`lmf`・`gwsc`

**6.2 MLO（局在軌道）が実用形に**: 既定 `mlo_method = 4`、パラメタは `mlo_delta`・`mlo_w`（eV、既定 2.0）だけ。`mlo_lm` で lm を原子ごとに指定
（旧 `Worb` は殻全体しか選べていなかった）。`Samples/MLOsamples` の 25 系を method 4 で回帰の試験に。FeMgO のスラブは空格子球 + ESM の 76 軌道の模型

**6.3 有限温度 QSGW（試験的）**: χ₀ 側 `t_tetrakbt`、Σ 側 `t_sigmaw`（2026-09-20 に `esmr`・`t_sigmakbt` を置き換え。`esmr` は同じ幅に換算して読む。
0.003 Ry → 262 K）。金属の QSGW が反復で荒れる正体は、小さい q の W_c(ω) の鋭いプラズモン極を Σ_c の実軸の極項が 1 点で拾うことだった
（LiTi₂O₄ で 1.8 eV、幅 0.1 eV）。`wcsmear = true`（既定）で中間の準位の均しを極の位置にも使うと、一発 GW の針（−3.4 eV）が消え、
9³ の反復の荒れが 48 → 13 meV。落とし穴: Σ_c の虚軸積分の Gaussian の正則化（PRB 76, 165106 式 57）は準位の均しそのもので、
極の項だけ核を変えると NiO 2³ の O 2s の対で 0.6 eV ずれた。虚軸の側も同じ Fermi–Dirac の核にした

**6.4 一発 GW と GPU**: `gw_lmfh <sname> -np N --gpu --mp`（`--mp` だけでは全ての積が TF32、精度が要るときは `--fp32`）。`hsfp0_gpu`（W の縮約の CUDA 化）

**6.5 直したバグ**: `lmf --quit=band` が LDA+U の `dmats` を NaN で壊していた／`zhev` は LAPACK の ier を先に見て原因を名指し、`zhev_tk2` の nev を n に丸める／
PROCAR の k 点の順（11 ランク以上で接尾辞の数値ソート）／`m_sxcf_sc` の実軸の極のビニングの範囲外アクセス、`readbandedge` の黙った代用／
`heftet` が絶縁体で `EFERMI_kbt` を書かなかった／MLO の `nskip` を k ごとに決めていて射影子が入れ替わり、バンドに折れが出た（全 k の最小に固定）／
`lmf` が前の実行の混合の履歴 `__mixm` を継いで、非磁性の解へ落ちた（Fe 2.13 → 0.02 μB。起動時に捨てる、`--keepmixm` で従来どおり）／
浅い半内殻の局在軌道（Zn 3d の `pz`）を MLO の模型関数に使えず、ZnO の MLO が 476 meV 外れた

**6.6 サンプルと試験**: `testecalj --all` は 26 ターゲット・64 件（2026-09-27）。`Samples/GetStarted/GaAs` が構造から QSGW のバンドまでの最短の道

## 7. ジョブを自動で回す仕組みの試作（`jobauto/`、2025-05）

いまの `ecalj_auto`（GW1500 の量産）の前の試作。kugui（ISSP）・ucgw を想定。考え方:

- 制御は手元で行う。計算機には ecalj を入れて PATH を通し、ssh で qsub と rsync ができればよい
- 遠くのディレクトリから必要なファイルだけを手元へ rsync して判断し、`qsub.0` を雛形に次のジョブを順に投げる。
  QSGW の反復ごとに作ると、物質の数 × 反復の数の `qsub.<id>` ができる
- 終わったかは、`qsub.0` が作るファイルで判断し、`qstat` は使わない
- kugui では `.bashrc` で `/etc/profile` を読む（ジョブの中でモジュールが使えるように）
- 例: `testSGA`（`init0/`）、`testSL3`（`init0/` に 170 ファイル）

## 8. 退役したスクリプト（`SRC/exec_legacy/`、2026-06-04 に退役、2026-10-01 に trash）

インストールされず、どこからも呼ばれていなかった。取り出すときは `git show c2d9df4aa:SRC/exec_legacy/<名前>`（表 1 の次のコミットの前）。

| 種類 | スクリプト |
| --- | --- |
| Wannier（MLWF）と相互作用 | `genMLWF2`・`genMLWF_vw`・`genMLWFhso`・`genMLWFmod`、`genuufpmt`、`jobmlw`、`FLEX_interaction.py`（cRPA による Hubbard の相互作用）、`RPAWannierTable.py`（2015、python2）、`Read_Sakakibara_W.py`、`Cal_W_off-site.py` |
| 誘電関数・磁気 | `epsPP0`・`epsPP0phispinsym`・`epsPP_lmfh`・`epsPP_lmfh_intra`・`epsPP_lmfhgamma`・`eps_lmfh`・`eps_lmfh_chipm`、`epsPP_magnon_chipm_mpi`、`chipm_mat`（2018）、`readeps.py`・`readepsNiSe.py` |
| GW のほかの計算 | `gwsigma`（Σ(ω)）、`gwgreen`、`gwphonon`、`job_sigm_fbz`（全 BZ の `sigm.<sname>.nk`、`--wsig_fbz` より精度がよい）、`seperate_SEComg*.py` |
| バンド・PDOS | `job_ham`・`job_ham_phispinsym`、`job_pdosc`・`job_pdosw`（`--ylmc` で複素 Ylm）、`jobmlo`、`ctrlgenPDOStest.py`（R は Å） |
| GWinput から TOML への移行の道具 | `tomlexpand.py`、`round_trip_check.sh`、`gwinput_xcheck.py`、`gwinput_unit_test.sh`、`check_defaults.py`、`disable_legacy_gwinput.py` |
| その他 | `Makeinstall`（CMake 以前のインストール）、`gwutil.py`、`a2vec.py`、`qqm`、`uutest` |

## 9. 古いサンプル（最上位の `MATERIALS/`、〜2026-05。2026-10-01 に `Samples/MATERIALS/` へ移した）

旧形式（`ctrl.<sname>`・`GWinput`）の入力と結果。どこからも参照されていなかった。trash ではなく `Samples/MATERIALS/` の下へ移した（user「いったん Samples の下へ」）。
近い内容のサンプルは `Samples` にある（`TestHomoDimerAtom` → `Samples/AtomDimer/N2`、`LaGaO3_relax` → `Samples/Relax/LaGaO3`、`Si_doping_sample` → `Samples/Doping/Si`、
`mass_fit_test` → `Samples/EffectiveMass`、`TETRAHEDRON_HomoGas` → `Samples/HomoGas`、`MLOsamples` → `Samples/MLOsamples`（旧形式の古い版）、
`cugase2_gwsc222`・`pdo_gwsc443`・`yh3fcc_gwsc666` → `Samples/TestInstall`、`GdNldau`・`erasldau` → 同じ系の LDA+U の試験 `TestInstall/gdn`・`eras`）。

- **原子の全エネルギーは一つに決まらない**（`Li_atom`）: (1) `lmfa` の球対称のスピン分極した原子の `etot`、(2) 胞に入れた球対称の原子の
  `lmf` の 1 回目の `ehf`、(3) 胞に入れた非球対称の原子の `lmf` の収束値。本来は (3) だが、胞と基底の打ち切りの誤差がある。
  (3) + ((1) − (2)) が補正したもの（(1) は球対称なら精度がよい）。一辺 6 Å では補正が 0.057 Ry と大きく、胞を大きくすると小さくなる
- `Fe_5`: Fe の Wannier と遮蔽した相互作用（`Coulomb_v.UP`・`Screening_W-v.UP`、25 MB）
- `NiMnSb_magnon`: NiMnSb のマグノン（MLWF 用と magnon 用の GWinput、`chipm_mat`、`memo_for_magnon`）
- `SiSigma`・`SiSigmaAny`: Σ(ω) の計算（`gwsigma`）と、任意の q（`<QPNT>` の AnyQ）
- `Si_HamMTO`・`Fe_HamMTO`・`FeMTOHAM`: MTO のハミルトニアン（`job_ham`）
- `NiSe_aftest`: 反強磁性の対称性の試験（`gwsc_sym`）。いまは `Samples/AFsymmetry/NiSe`
- `TESTsamples`（157 ファイル、2012 年頃の ctrl の例）、`SYMLsamples`（`syml` の例 22）、`CTRLsample`（Cu・Fe・Si の ctrl）、
  `EPS_Ag`、`LiDOS_Discrete`、`H_atom`、`NiO`、`Materials.ctrls.database`・`job_materials.py`・`job_vbm`・`sortvbm.py`（多くの物質を回す道具）

---

## 表 1. 片付けたもの（trash に移したもの）

最上位の `MATERIALS/` は trash ではなく `Samples/MATERIALS/` へ移した（§9）。

| 日 | もの | 元の場所 | 外す前のコミット | 過去ログの節 |
| --- | --- | --- | --- | --- |
| 2026-10-01 | Materials Project の API キーを含む設定の写し | `ecalj_auto/OUTPUT/*/config.ini` の `apikey` 行 | `290397b34` の前 | ecalj_auto/README.md |
| 2026-10-01 | 古い写し | `SRC/BK`、`SRC/execgfortran`、`SRC/execAHC` | `116254e1d` の前（`4f9332d98`） | — |
| 2026-10-01 | GW1500 の古いスクリプト・表 | `ecalj_auto/run_gw1500_addrun*.sh`、`jobgw1500.sh`、`gw1500_recheck_*.sh`、`qpu_change.py`、`gw1500_rerun_20260930.tsv` | `116254e1d` の前 | ecalj_auto/GW1500_status.md |
| 2026-10-01 | ビルドの生成物 | `SRC/exec/build`・`build_gf14`、`SRC/exec_gfortran`・`exec_gfortran-14`・`exec_ifx` | `3f0771f2f` の前 | §1.5 |
| 2026-10-01 | 打ち込み用・古いスクリプト | `g`・`gg`・`i`・`ii`・`n`、`build_nvfortran.sh`、`gpu_wait_run.sh`、`jobinstall_kugui.sh`・`jobinstall_ohtaka.sh` | `3f0771f2f` | §1.1〜1.3、§2.1 |
| 2026-10-01 | 古いビルドとソースの控え | `SRC/exec/BK` | `3f0771f2f` | §1.4 |
| 2026-10-01 | リファクタの調査と回帰の控え | `.refactor_notes/` | `3f0771f2f` | §3.1 |
| 2026-10-01 | サンプルの古い試行と控え | 表の下の注 | `3f0771f2f` | §4.1、§4.2 |
| 2026-10-01 | hgw の統合の記録 | `PHASE1B_REFACTOR.md` | `c2d9df4aa` | §3.2 |
| 2026-10-01 | ecaljdoc の古い下書き | `ecaljdoc_drafts/` | `c2d9df4aa` | §4.3 |
| 2026-10-01 | 一発 GW と有限温度の手引き | `FiniteT_and_QPE_HOWTO.md` | `c2d9df4aa` | §5 |
| 2026-10-01 | 2026-06〜09 の更新の要約 | `HIGHLIGHTS_2026-06_09.md` | `c2d9df4aa` | §6 |
| 2026-10-01 | ジョブの自動実行の試作 | `jobauto/` | `c2d9df4aa` | §7 |
| 2026-10-01 | 退役したスクリプト | `SRC/exec_legacy/` | `c2d9df4aa` | §8 |

注: サンプルの古い試行と控えは、MLOsamples の `test*`・`temp`・`*.bk`・`*.tmp`（§4.1）、`Samples/TestInstall/TESTunused`、
`Samples/TestInstall/eras/occnum.eras.bk`、`TOOLS/FparserTools/f_calltree.py.bk*`、`TOOLS/SrcFragments/f_calltree.py.bk*`・`ANALYZEnotusednow/analyze_temp~`、
`ecalj_auto/INPUT/testSGA/joblist.bk`。
