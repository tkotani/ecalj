# `lmf` / `lmfa` / `lmchk`（とその仲間）のコマンドラインのオプション

## コマンドラインの引数の形（最小限）

Fortran のプログラムは、次の形の引数を受け付けます。それ以外は、起動したところで一行のヒントを出して止まります。
引数の解析はすべて
[`SRC/subroutines/m_cmdopt_registry.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/m_cmdopt_registry.f90)
（オプションの表、打ち間違いと廃止した書き方の検査）と
[`SRC/subroutines/m_toml_override.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/m_toml_override.f90)
（TOML の上書き）にあります。

表 1. 引数の形

| 形 | 意味 | 行き先 |
|---|---|---|
| `<sname>` | 位置引数。ctrlg のファイルの拡張子 | `m_ext::sname` |
| `ctrlg.<sname>.toml` | 位置引数。TOML のファイル名そのもの（TAB で補完したときなど） | `m_ext::sname` |
| `--ctrlg:<ドットで区切ったパス>=<値>` | 実行時の TOML の上書き（メモリの中で。ファイルは変えない） | `m_toml_override` |
| `--<flag>` | 登録された cmdopt0。あるかないかの真偽 | `m_cmdopt_registry` の `c0_<flag>` |
| `--<flag>=<value>` | 登録された cmdopt2。型のある値 | `m_cmdopt_registry` の `c2_<flag>` |

二つの位置引数の形は同じ意味です。`lmf nio` も `lmf ctrlg.nio.toml` も `sname='nio'` になります。
`m_ext_init` はその中間の形を受け付けません。`ctrlg.nio`（`.toml` が無い）と `ctrl.nio`（古いテキストの形式）は、
二つの正しい形（後者には `Legacy2toml.py`）を示すヒントを出して止まります。

起動のときに、引数を一度だけ検査します（`m_cmdopt_registry::validate_arglist`、`m_args::m_setargs` から呼ぶ）。次の場合に止まります。
- cmdopt0 の表に無い `--<word>` / `-<word>`
- 頭の部分が cmdopt2 の表に無い `--<word>=<value>`
- 下の廃止した形のどれか（形ごとに移行のヒントを出す）

新しいオプションを足すのは、`m_cmdopt_registry.f90` にそれぞれ一行です。`c0_<name>`（または `c2_<name>`）を宣言し、
該当する `load_*_registry` に `call set0('--<name>', c0_<name>, narg, arglist)`（または `get2('--<name>', outs, narg, arglist)`）を足します。
打ち間違いを見つける表は、その副産物として作られます。別に「知っているフラグの一覧」を持って合わせる必要はありません。

このページは、登録されたオプションの一覧です。手で書いていて、
`grep -E "call set0\('--|get2\('--" SRC/subroutines/m_cmdopt_registry.f90`
（または `lmf --listcmdopt`）で確かめています。

## 表の読み方

- **フラグ**は cmdopt0（値なし）で、そのまま書きます。例 `--writeham`。
- **キー=値**は cmdopt2 で、必ず `=` を付けます。例 `--jobgw=1`。知らない `--<word>=<value>` は打ち間違いのヒントを出して止まります。
- **TOML の上書き**は別の名前空間です: `--ctrlg:<ドットで区切ったパス>=<値>`（`m_toml_override.f90` が扱う）。この表には載せません。[下](#run-time-ctrlg-overrides) を参照。
- ほとんどのフラグは `lmf`（と、GW の準備をする `lmf --jobgw=…`）で使います。`lmfa` / `lmchk` が使うものは、そう書いてあります。

## 廃止した形（移行のヒントを出して止まる）

以前は TOML のキーを直接上書きする `cmdopt0` / `cmdopt2` の近道でした。いまはどれも、コマンドラインにあると、
正しい形 `--ctrlg:<path>=<value>` を示す一行を出して止まります。

表 2. 廃止した形と、代わりの書き方

| 廃止した形 | 代わり |
|---|---|
| `--pr=N` / `--pr N` / `-pr N` | `--ctrlg:verbose=N` |
| `--time=N[,M]` | `--ctrlg:time=[N,M]` |
| `--phispinsym`（コマンドラインのフラグ） | `--ctrlg:ham.phispinsym=true` |
| `-v<name>=<value>`（古い `%const`） | `--ctrlg:<dotted.path>=<value>` |
| `-v[<path>]=<value>`（古い角かっこの形） | `--ctrlg:<dotted.path>=<value>` |
| `--[<path>]=<value>`（移行の途中の形） | `--ctrlg:<dotted.path>=<value>` |
| `--<a.b.c>=<value>`（移行の途中の形） | `--ctrlg:<a.b.c>=<value>` |
| `--toml.<path>=<value>`（移行の途中の形） | `--ctrlg:<path>=<value>` |
| `-emin=` / `-emax=` / `-ndos=`（ダッシュ一つ、古い形） | `--emin=` / `--emax=` / `--ndos=` |
| `-EfermiShifteV=`（ダッシュ一つ、古い形） | `--EfermiShifteV=` |
| `--terse` / `-terse`（lmchk の重なりの表） | なくした。`--ctrlg:verbose=N`（N≤40）を使う |
| `--no-iactiv`（古い lmfa のフラグ） | なくした。以前は黙って無視していたが、いまは止まる |

## 共通（lmf / lmfa / lmchk）

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--help` | flag | 使い方とこのページへのリンクを出して終わる | `m_ext.f90` |
| `--debug` | flag | 全体のデバッグ出力 | 多数 |
| `--show_time` | flag | `hambl` の時間を出す | `hambl.f90` |
| `--fullstdo` | flag | ランク 0 だけでなく、すべてのランクが出力する。GW のプログラムでは、各ランクが自分の `stdout.<rank>.<prog>` を書く（`gwsc` はそれを `STDOUT/` に移す）。無いときは、ランク 0 が実行のログに書き、ほかのランクは /dev/null に書く（2026-09-27） | `m_mpi.f90` |

## 自己無撞着の計算・電子状態（lmf の本体）

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--etot` | flag | 全エネルギーの診断を出す | `m_lmfinit.f90` |
| `--vbmonly` | flag | 真空に対するおおよその VBM/CBM を出して終わる | `main_lmf.f90` |
| `--getq` | flag | チャネルごとの占有 Q を出す | `main_lmf.f90` |
| `--v0fix` | flag | 自己無撞着の計算の間、球対称なポテンシャル V₀ を固定する | `m_lmfinit.f90` |
| `--keepmixm` | flag | 前の `lmf` の混合の履歴 `__mixm.<sname>` を残す（無いと `lmf` は始めに捨てる） | `m_lmfp.f90` |
| `--noinv` | flag | BZ のサンプリングで反転対称を使わない | `m_lmfinit.f90` |
| `--nosym` | flag | 空間の対称性を使わない（対称性なしの計算にする） | `m_lmfinit.f90` |
| `--nosymdm` | flag | LDA+U の密度行列に対称性を使わない | いろいろ |
| `--shorten` | flag | 格子ベクトルを短くする（`lmchk` の動作） | `lmaux.f90` |
| `--slat` | flag | Ewald の短くした格子を出す（`lmchk`） | `lmaux.f90` |
| `--afsym` | flag | 反強磁性（Shubnikov）の対称性を使う。BZ のサンプリングが半分になる | `m_qplist.f90` |
| `--onesp` | flag | `nsp=2` でも `nspxx=1` にする | `m_qplist.f90` |
| `--avoidgamma` | flag | q のサンプリングで Γ 点を飛ばす | いろいろ |

## 途中で止める点（デバッグ、段階を分けた実行）

| オプション | 働き | 使う場所 |
|---|---|---|
| `--quit=show` | `--show` の表示の後で止まる | `m_lmfinit.f90` |
| `--quit=ham` | `m_hamindex_init` の後で止まる | `main_lmf.f90` |
| `--quit=mkpot` | `m_mkpot`（ポテンシャルを作る）の後で止まる | `m_bndfp.f90` |
| `--quit=dmat` | 最初に密度行列を作った後で止まる | `main_lmf.f90` |
| `--quit=band` | バンドのループの後で止まる（密度は更新しない） | `m_lmfinit.f90` |
| `--quitecore` | コアの扱いの後で止まる | `m_bndfp.f90` |

## バンド・状態密度・フェルミ面

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--band` | flag | `syml.<sname>` に沿って固有値を計算する（バンドプロット） | `m_lmfinit.f90` |
| `--fullmesh` | flag | BZ の全メッシュ。フェルミ面、PROCAR、BoltzTraP で使う | `m_qplist.f90` |
| `--fermisurface` | flag | フェルミ面のモード（XCrySDen の `.xsf`） | `m_bndfp.f90` |
| `--mkprocar` | flag | PROCAR（原子・軌道へ射影した重み）を書く | `m_bndfp.f90` |
| `--writepdos` | flag | `pdos.*` のファイルを書く | `main_lmf.f90` |
| `--dos` / `--pdos` / `--tdos` | flag | 状態密度 / 部分状態密度 / 全状態密度 | `m_lmfinit.f90` |
| `--tdostetf` | flag | `--tdos` で四面体法を使わない（代わりに Methfessel–Paxton） | `m_bndfp.f90` |
| `--wdsawada` | flag | Sawada の形の簡単な状態密度を出して終わる | `main_lmf.f90` |
| `--eigen-at-k` | flag | 与えた k での固有値を書く | `m_qplist.f90` |
| `--writeeigen` | flag | `--eigen-at-k` と同じ（別名） | `m_bndfp.f90` |
| `--allband` | flag | すべての k ですべてのバンドを出す | `m_bandcal.f90` |
| `--boltztrap` | flag | 全メッシュで `BoltzTraP` の入力を書く | `m_bndfp.f90` |
| `--cls` | flag | コアの準位の分光の重み（`m_clsmode.f90`） | `m_bndfp.f90` |
| `--corehole` | flag | コアホールの計算 | いろいろ |
| `--efermi=<file>` | kv | `lmf` が書き（バンドプロットのモードでは読み）、`mlo` が読むフェルミ準位のファイルを、`efermi.lmf` の代わりにこれにする。`job_band`、`job_dos`、`job_fermisurface` は `efermi.lmf.job_band` などを、`job_mlo_soc` は `efermi_soc` を渡す。こうして `efermi.lmf` には最後の自己無撞着な `lmf` の値が残る（2026-09-26） | `m_cmdopt_registry.f90`（`efermi_file`） |

## GW の準備とその仲間（`lmf --jobgw=…`、`gwsc` の中）

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--jobgw=0` | kv | `lmf` が索引のファイル `__HAMindex0`（`qg4gw` と `gwinit` のため）を書いて止まる | `main_lmf.f90` |
| `--jobgw=1` | kv | `lmf` が GW の入力を作る（古い `lmfgw-MPIK` の代わり）。`hgw --jobgw=1` は `gwsc` の GW の段 | `main_lmf.f90`、`main_hgw.f90` |
| `--writeham` | flag | 後の道具（mlo / writeham）のために `HamiltonianPMT.*` を書く | `main_lmf.f90` |
| `--writev0` | flag | V₀（球対称部分）を動径メッシュで書く | `locpot.f90` |
| `--writesene` | flag | 自己エネルギーのファイル `__sene.*` を書く（デバッグ） | `m_bndfp.f90` |
| `--wsig_fbz` | flag | BZ 全体に並べ直した `sigm` を書く | `m_bndfp.f90` |
| `--use_sigm_fbz` | flag | ふつうの `sigm.*` の代わりに `sigm_fbz` を読む | `m_bndfp.f90` |
| `--readQforGW` | flag | `qg4gw`: `[gw] QforGW` の q ベクトルを使う（one-shot の GW。`gw_lmfh` が渡す） | `m_q0p.f90` |
| `--n1n2n3eps` | flag | `qg4gw`/`mkqg` で、ε/χ⁺⁻ のために `n1n2n3eps` のメッシュを使う | `mkqg.f90` |
| `--job=N` | kv | `hbasfp0`、`hvccfp0`、`hsfp0_sc`、`hsfp0`、`hx0fp0`、`hqpe`、`qg4gw`、`heftet` の段を選ぶ。範囲はプログラムごと。`gwsc` が使う値は [gwsc](./gwsc) の画面の出力にある | 多数 |
| `--ntqxx` | flag | 各 q の Σ のバンドの数を、最初の反復の値に固定する（ファイル `NTQXX`。[gwsc](./gwsc)） | `m_itq.f90` |
| `--nk=` | kv | `mlo_magnon` の k の並列のグループの大きさ（ランクの数から決まる既定を上書き） | `main_mlo_magnon.f90` |
| `--cutuu=`、`--nww=`、`--nb=`、`--Wtype=` | — | 廃止（2026-10-02、AHC・Wannier のマグノンを外した）。与えると止まる | `m_cmdopt_registry.f90` |
| `--dumpW` | flag | `hgw`: メモリに持つ W−v を、解析のために `__WVR.<iq>` / `__WVI.<iq>` に書く（`gwsc --dumpW`） | `main_hgw.f90` |
| `--WVR2ptRaxis` | flag | 実軸での Σc の極の項: 極での W_c を 3 点でなく 2 点の線形内挿で求める（`gwsc --WVR2ptRaxis`） | `m_sxcf_sc.f90` |
| `--wcsmear` | flag | `[gw] wcsmear = true` と同じ（キーの方を使う。既定で true） | `m_sxcf_sc.f90` |
| `--sigma_tf32` | flag | `--use_fp32` の MP の GPU のプログラムで、Σc の積だけを仮数 10 ビットの入力で行う（`gwsc --prec=tf32` が渡す）。[GPU 版](./ecaljgpu) を参照 | `m_sxcf_sc.f90`、`m_linalg_policy.f90` |
| `--tetwt_write` | flag | `hgw`: すべての q の四面体の重みを `__TETWT.<iq>.<isp>` に書くだけで止まる（GPU、W、Σ なし）。`gwsc --gpu` が `hgw` の横の CPU のコアで走らせ、`hgw` はそれを読む（ビットまで同じ。無いか合わないファイルは `hgw` が計算する） | `main_hgw.f90`、`x0kf_v4h.f90` |
| `--qibzonly` | flag | 既約 BZ の q 点だけを計算する（補助の q₀ を飛ばす） | いろいろ |
| `--geteta` | flag | ε(q,ω) の η の減衰の診断 | いろいろ |
| `--interbandonly` / `--intrabandonly` | flag | ε でバンド間・バンド内の遷移だけにする | いろいろ |
| `--testso` | flag | SO の確認の経路 | いろいろ |
| `--tetwtk` | flag | χ0: 四面体の重みを、すべての k をまとめてでなく、k のループの中で k ごとに計算する（メモリが少ない） | `x0kf_v4h.f90` |
| `--tetraw` | flag | `lmf` の PDOS のモード: 四面体のデータ `tetradata.dat` と k の一覧 `qplistf.dat` を書く | `m_procar.f90` |

## MLO（Wannier の代わり）（`job_mlo`、`job_mloW`、`job_mlo_soc`）

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--mlo` | flag | MLO のモード。`HamRsMLO` があるとき、`lmf --jobgw=1 --mlo` は `__cmlo.data`/`.info` を書く（段 a'）。`mlo --mlo` は `HamRsMLO`（MLO の索引）を作り、`__HamiltonianGW.info` があれば（job_mloW、job_mlo_magnon）そこから `__cmlo` も作る | `sugw.f90`、`m_HamPMT.f90` |
| `--mlog` | flag | MLO と GW をまとめた出力 | いろいろ |
| `--mlo_diagnorm`、`--mlo_feb4`、`--mlo_ortho`、`--mlo_orthonorm`、`--gs` | — | 廃止（2026-10-02）。与えると止まる。MLO は実空間で 2 乗積分が 1 になるよう、軌道ごとに k によらない定数で規格化する（[mlo](./mlo) §6）。k ごとの規格化と直交化は、軌道の形と内挿したバンドを変えるので使わない | `m_cmdopt_registry.f90` |
| `--mlofreeze` | flag | MLO-gwsc（開発中）: 既存の `HamRsMLO` の MLO の索引を使い、`QMLO_SigRs` だけを書き直す。`gwsc --mlo` が反復の中の `mlo` に渡す。[MLO-gwsc](./mlo_gwsc) を参照 | `m_HamPMT.f90` |
| `--dwnb=mlo` | kv | `huumat_MPI`: ⟨u_k\|u_k+b⟩ の基底を MLO にする（`job_mlo_magnon`、`mlo_spread.py`）。`--dwnb=wan`（Wannier 関数）は 2026-10-02 に外した | `main_huumat.f90` |
| `--ahc`、`--mloahc`、`--AHCMAT`、`--UUMAT`、`--cmlo`、`--q2q1test` | — | 廃止（2026-10-02、AHC・lmfham2・試験用）。与えると止まる | `m_cmdopt_registry.f90` |
| `--writedw` | flag | 既定でオン: `m_procar` で `dw` の重みを書く | `m_procar.f90` |
| `--nowritedw` | flag | `dw` を書かない | `m_procar.f90` |
| `--socmatrix` | flag | SOC の行列（`__HamiltonianPMTsoc`）を書く | `m_HamPMT.f90` |

## 密度・ポテンシャルの入出力

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--density` | flag | なめらかな密度を `xsf` などのファイルに書く（`m_bndfp.f90`） | `m_bndfp.f90` |
| `--wrhomt` | flag | MT の中の本当の密度 `rhoMT.<ib>` を書く（デバッグ） | `locpot.f90` |
| `--wpotmt` | flag | MT の中の本当のポテンシャル `vtrue.<ib>` を書く | `locpot.f90` |
| `--vesatom` | flag | `lmfa` の静電ポテンシャルを出す | `freeat.f90` |
| `--vesdat` | flag | ランクごとの静電のデータを書く | `locpot.f90` |
| `--espot` | flag | 静電ポテンシャルの計算 | `m_bndfp.f90` |
| `--estaticall` | flag | 静電ポテンシャル全体を出す | `m_bndfp.f90` |
| `--eszero` | flag | Madelung の補正を 0 にする（診断） | いろいろ |
| `--novxc` | flag | 交換相関の寄与を飛ばす（デバッグ） | `locpot.f90` |

## 構造・対称性（`lmchk` の動作）

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--getwsr` | flag | 接する MT 半径を求める。`ctrlgenToml.py` が使う | `lmaux.f90` |
| `--kchk` | flag | k のメッシュが空間群と合うかを確かめる | `bzmesh.f90` |
| `--ylmc` | flag | `Ylm` の解析を出す | いろいろ |

## SOC

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--socmatrix` | flag | （MLO の節を参照） | `m_HamPMT.f90` |
| `--skiphammsoc` | flag | `hamm` に SOC の項を足さない（摂動の経路） | `m_bandcal.f90` |
| `--ctrlg:ham.phispinsym=true` | TOML | スピンで平均した動径関数（以前の `--phispinsym`）。TOML の上書きか、`ctrlg.<sname>.toml` の `[ham]` で設定する | `m_toml_override.f90` |

## GPU・線形代数の調整

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--gpu` | flag | GPU のコードを使う（`InstallAll.py` は `/tmp/gpu.lock` を取り、`gwsc` は `/tmp/ecalj_res/` に GPU ごとのロックを取る。[GPU 版 § 複数のジョブ](./ecaljgpu#gpu-multijob)） | `InstallAll.py`、`m_*.f90` |
| `--linalg=<row>,<row>,…` | kv | GPU の行列演算の表の行を上書きする。例 `--linalg=fp32.cgemm.large=realsgemm,fp32.epsinv.large=mixed1`。表と方法は [GPU 版 § 自動選択](./ecaljgpu#gpu-precision) を参照 | `m_linalg_policy.f90` |
| `--use_gemmul8` | flag | 大きな GPU の行列の積を GEMMul8（INT8 のテンソルコアでの Ozaki-II）で行う。`InstallAll.py --gemmul8` が要る。[GPU 版](./ecaljgpu) を参照 | `m_blas.f90`、`m_gemmul8.f90` |
| `--modifiedGS` | flag | 修正 Gram–Schmidt の直交化 | いろいろ |
| `--diag=default` / `--diag=tridiag` / `--diag=chefsi` | kv | 対角化の方法を選ぶ（LAPACK の既定、三重対角、ChebFSI） | いろいろ |

## 後処理の道具（`hmlo_ovlppair`、eps の道具）

| オプション | 型 | 働き | 使う場所 |
|---|---|---|---|
| `--sp1=` / `--sp2=` | kv | `hmlo_ovlppair` のスピンの選択 | `hmlo_ovlppair.f90` |
| `--emin=` / `--emax=` / `--ndos=` | kv | 状態密度のエネルギーの窓と分解能 | いろいろ |
| `--EfermiShifteV=` | kv | 状態密度やスペクトルの前に `EF` をずらす（eV） | いろいろ |

## 診断・デバッグだけのフラグ

開発のためのもので、ふつうの計算では役に立ちそうにありません。振る舞いが変わらない保証もありません。

`--debugbndfp`、`--debugpwmat`、`--debugsugw`、`--debugzmel`、`--zmel0`、
`--normcheck`、`--showdmat`、`--x0test`、`--cvK:`、
`--skip_qvalcheck`、`--skipbstruxinit`、`--skipCPHI`、
`--skiphammsoc`、`--skip1d` / `--skip2nd` / `--skip2ndd` / `--skip2ndp`
/ `--skip2nds`、`--skipd` / `--skipf` / `--skiplo`（どれもいまは `m_HamPMT.f90` の使われないコード）、`--wanatom`。

`--skipGS` はデバッグ用ではありません。`lmf --jobgw=1`（`sugw.f90`）の Gram-Schmidt の規格直交化を飛ばして速くします。[gwsc](./gwsc) を参照。

## 実行時の TOML の上書き（`--ctrlg:<path>=<value>`） {#run-time-ctrlg-overrides}

上書きの書き方は `--ctrlg:<ドットで区切ったパス>=<値>` だけです。
値は TOML の型で読みます（真偽は小文字の `true`/`false`、文字列は引用符で囲む、配列は `[...]`）。
裸の識別子には自動で引用符を付けるので、シェルが引用符を取ってしまっても困りません。

```bash
lmf si --ctrlg:bz.nkabc=[8,8,8] --ctrlg:ham.scaledsigma=0.8
lmf si --ctrlg:ham.so=1 --ctrlg:ham.nspin=2          # SOC as perturbation
lmf si --ctrlg:ham.phispinsym=true                   # bool lowercase
lmf si --ctrlg:verbose=50 --ctrlg:time=[5,5]         # top-level keys
lmf si --ctrlg:spec.1.r=2.5                          # override SPEC[1].r
lmf si --ctrlg:symgrp=find                           # bare identifier auto-quoted
```

パスは TOML の構造と同じです（`section.key`。表の配列では `section.<idx>.key` で、番号は 1 から）。
使った上書きはランク 0 がログに書きます。ファイルに無いキーは書き足されます
（その節の見出しの下に。最上位のキーはファイルの先頭に。節が無ければ新しい `[section]` を作って）。そのこともログに出ます。
2026-09-30 までは、そのような上書きは飛ばされ、既定値のまま計算が進みました。
綴りの違うキーも書き足され、どのプログラムも読まないので、ログの行で確かめてください。
`--ctrlg:` の説明はこの節が正本です。古い `ctrl.<sname>` の `%const` の定数は、変換のときに `ctrlg.<sname>.toml` に書き込まれる（[入力ファイル](./toml_migration)）ので、TOML のパスで上書きします。

## このページが古くなったとき

ソースにあってこのページに無いフラグ（またはその逆）があれば、ソースを信じてください。次で調べ直せます。

```bash
lmf --listcmdopt | sort        # the registered names (cmdopt2 with a trailing =)
```

そのうえで手で合わせます。このページを自動では作りません。分類と説明の文に価値があるからです。アルファベット順の一覧なら `grep` で得られます。

## TAB による補完（bash）

`InstallAll.py` は、cmdopt の登録と TOML の上書きの書き方の両方を知っている bash の補完のスクリプトを入れます。インストールの後は:

- `lmf <TAB>`（または `lmfa`、`lmchk`）→ いまのディレクトリの `ctrlg.<sname>.toml` を、位置引数の `<sname>` として出す。
- `lmf nio --w<TAB>` → `--w` で始まる cmdopt0/cmdopt2 の名前（`--writedw`、`--writeham`、`--writepdos`、...）。
- `lmf nio --ctrlg:<TAB>`（または `--ctrlg:bz.` のような途中まで）→ いまのディレクトリの `ctrlg.<sname>.toml` に**実際にある**ドットで区切ったキー。
  `[[spec]]` / `[[site]]` の配列は、自動で `spec.1.r`、`spec.2.r`、... に展開する。
- `--ctrlg:` と cmdopt2 の後ろには空白を付けない（`--ct<TAB>` は `--ctrlg:` になり、カーソルはコロンの後ろに残る。続けて `bz.nkabc=...` と書ける）。

補完のデータの出どころは二つで、新しさが違います（表 3）。

表 3. 補完のデータの出どころ

| 補完するもの | 出どころ | 更新されるとき |
|---|---|---|
| `--ctrlg:<dotted.path>` | いまのディレクトリの `ctrlg.<sname>.toml` を、TAB のたびに `python3`（tomllib）で読む | いつもディスクの TOML のとおり |
| `--<flag>` / `--<flag>=` | `BINDIR/ecalj_cmdopts.list`（`InstallAll.py` のときに `lmf --listcmdopt` が書く） | `./InstallAll.py` を走らせるたび |

対応は解析の仕組みと同じです。

- `ctrlg.*.toml` のキーは、そのファイルに本当にあるときだけ補完に出ます。そのため、TAB で補完した `--ctrlg:` のパスは綴りの違うキーになりません（ファイルに無いキーは書き足される。2026-09-30 から）。
- `--<flag>` は、`m_cmdopt_registry::load_cmdopt0_registry` / `load_cmdopt2_registry` に登録されているときだけ出ます。そのため、TAB で補完したフラグが `validate_arglist` の「知らないオプション」で止まることはありません。

`~/bin/ecalj_complete.bash` と `~/bin/ecalj_cmdopts.list` は、どちらも `InstallAll.py` が書きます。
`~/.bashrc` には、`[ -f .../ecalj_complete.bash ] && source ...` の一行を一度だけ足します。
足したくないときは `InstallAll.py --no-bashrc` とします。あとで切るときは、`~/.bashrc` の `# >>> ecalj bash completion ... <<<` の部分を消します。

`lmf --listcmdopt`（補完の仕組みと一緒に足した）は、登録されたすべての cmdopt0 / cmdopt2 の名前を、一行に一つずつ標準出力に出します。
cmdopt2 には後ろに `=` を付けるので、補完のスクリプトはただのフラグと区別できます。

```text
$ lmf --listcmdopt | head -3
--afsym
--allband
--avoidgamma
$ lmf --listcmdopt | grep '=$' | head -3
--jobgw=
--job=
--nk=
```

インストールのときに登録の中身を知る必要のある、外のスクリプトのデータも、同じ `--listcmdopt` で作れます。
