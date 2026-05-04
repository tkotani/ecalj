# ecalj Claude Development Guide

ecalj を Claude (LLM) で開発する際のガイドライン。
コーディング規約、ビルドの罠、GPU 開発の教訓、アーキテクチャ方針をまとめる。

## コーディング規約: Fortran Singleton Pattern (kotani 設計)

ecalj は **Fortran module を singleton class として使う** 設計パターンを採用している。
OOP の GoF singleton pattern を Fortran module の言語特性で自然に実現したもので、
科学計算コードにおいて `type` (derived type) を避けつつ、状態管理・データフロー明示・
GPU (OpenACC) 互換性を同時に達成する実践的手法。

### なぜ singleton か

Fortran の module は言語仕様上、プログラム内で唯一のインスタンスを持つ。
これは OOP の singleton と同じ性質であり、科学計算の多くのサブシステム
(ハミルトニアン、波動関数ストレージ、自己エネルギー計算器など) は
複数インスタンスが不要なため、module = singleton が自然に適合する。

C++ や Python では singleton pattern を「わざわざ」実装する必要があるが、
Fortran では module 自体がそれを提供する。ecalj はこの言語特性を
意識的・体系的に活用している。

### 設計原則
- **1 module = 1 責務 = 1 状態セット**: 複数インスタンス不要なら `type` を被せない
- **状態は module 変数で公開**: `protected, public` で読み取り公開、書き換えは module 内 subroutine 経由のみ
- **import は明示**: 必ず `use m_foo, only: bar, baz` で必要な symbol だけ取る
- **subroutine 規約**: 入力=引数、出力=module 変数で返す
- **`type` (derived type) は積極的に避ける**: Fortran の type には以下の罠がある
  - `intent(out)` で allocatable components が暗黙に deallocate される
  - `present(state%alloc_member)` が nvfortran + OpenACC で silent fail
  - 合成 type の copy semantics が不明瞭
  - `use only` でデータフローが見えない (type の中身は呼出側から不透明)

### 階層化設計: module が module を use する DAG

singleton module 間のデータ受け渡しは `use` の階層で表現する。
各 module は自分より下位の module だけを `use` し、循環依存を禁止する。
これにより依存関係が DAG (有向非巡回グラフ) になり、データの流れが一方向に追える。

```
m_struct_from_lmf (構造: plat, alat, pos, z)
    ↑ use
m_gw_product_basis (積基底: ndima, nlnx)
    ↑ use
m_sxcf_sc (自己エネルギー計算)
    ↑ use
main_hrcxq / hgw_combined (メインプログラム)
```

**データフローの規約:**
- **大量データの入力**: `use m_foo, only: bar` で下位 module の `protected` 変数を読む。subroutine 引数では渡さない
- **subroutine 引数**: フェルミエネルギー `ef`、q 点 `qp`、スピン `isp`、反復回数 `niter` など、**指示的・局所的な値に限る**。配列データを引数で渡し回さない
- **出力**: subroutine が計算した結果は module 変数に格納し `protected` で公開。caller は `use only` で読む

この規約により:
- `type(state)` を引数で渡し回す必要がない (type を避けられる)
- どの module がどのデータを提供し、誰が消費するかが `use only` で完全に見える
- `grep 'use m_foo' SRC/subroutines/*.f90` で依存関係を即座にトレースできる
- module の差し替え (例: FILE → MEMORY backend) が caller に影響しない

### なぜ Fortran OOP (type-bound procedure, class) を使わないか

Fortran 2003 で `class`, `type-bound procedure`, `select type`, inheritance が導入されたが、
これらは C++/Java のために設計された OOP 概念を配列計算特化言語に後付けしたもので、
実用上の問題が多い:

- `select type` による polymorphism は冗長で可読性が低い
- `type-bound procedure` は module procedure への syntax sugar に過ぎない
- inheritance は Fortran のメモリモデル (allocatable component) と相性が悪い
- コンパイラのバグが多い (特に nvfortran + OpenACC + derived type の組み合わせ)

ecalj の singleton + 階層 use は、Fortran の元々の強み (module, allocatable, 配列操作) だけで
設計を完結させる **言語に逆らわない** アプローチ。

### singleton の利点

1. **OpenACC/GPU との相性**: device 変数の管理が単純。module 変数は全 subroutine から直接アクセスでき、`!$acc` ディレクティブで素直に扱える。type のメンバを device に置くと nvfortran の制約に次々ぶつかる
2. **データフローが grep で追える**: `use m_foo, only: bar` を grep すれば誰が何を使っているか一目瞭然。type だと `state%bar` を追う必要があり、複数の type が絡むとトレースが困難
3. **re-entrant 化が楽**: module 変数を dealloc-realloc すればいい。type のネストした copy semantics を考えなくていい
4. **backend 切替が caller 透過**: module 内部で FILE / MEMORY_3D 等の実装を切り替えても、caller は同じ module 変数を `use` するだけ
5. **LLM にも優しい**: module のヘッダ (変数宣言部) を見れば状態が全部わかる。type の定義を辿って中身を理解する必要がない

### 核心的利点: 変数宣言の唯一性

従来の subroutine 引数渡し方式では、同じ変数を caller と callee の両方で宣言する
必要がある (二重記述):

```fortran
! 呼ぶ側
call calc_sigma(ndim, nqibz, nspinmx, qibz, wk, ekc, ...)

! 呼ばれる側 — 同じ変数を再宣言、サイズも再指定
subroutine calc_sigma(ndim, nqibz, nspinmx, qibz, wk, ekc, ...)
  integer, intent(in) :: ndim, nqibz, nspinmx
  real(8), intent(in) :: qibz(3,nqibz), wk(nqibz), ekc(ndim)
```

引数が 10〜20 個になると caller と callee の対応がずれてバグになる。
配列サイズの整合性も人間が保証しなければならない。二重記述は百害あって一利なし。

singleton + `use only` なら変数宣言は module 内の一箇所だけ:

```fortran
! 呼ばれる側 — use で取るだけ、再宣言不要
subroutine calc_sigma(ef, qp, isp)  ! 指示値のみ引数
  use m_read_bzdata, only: qibz, wk, nqibz
  use m_struct_from_lmf, only: nspinmx
```

- 変数の生成元が唯一 — 二重記述ゼロ
- サイズ不整合が原理的に起きない
- 引数リストが短い (指示値のみ) ので caller-callee の対応ミスが起きない

**深い呼び出し階層での追加データ問題**: 5段の call 階層の末端で `natom` が1つ必要になった場合:
- **引数渡し**: 途中の全 subroutine (5箇所) の引数リストに `natom` を追加。大量の変更、リグレッションリスク
- **type**: `state%natom` として渡すと、データが本来の出自 (生成元 module) から切り離されコピーが増える。実際に起きるのは、たまたま 5 段目まで通っている type に `nbas` を便乗で放り込んでしまうこと。やがて `crystal_state%nbas`、`gw_state%nbas`、`basis_state%nbas` と複数の type に同じ値のコピーが散在し、どれが正でいつ同期されるか誰にもわからなくなる — **データの唯一性が簡単に消失する**
- **singleton use**: `use m_struct_from_lmf, only: natom` を末端の1箇所に追加するだけ。途中の subroutine は変更不要、データの出自は `m_struct_from_lmf` と明確

### 弱点と対策

- **subroutine の出力がシグネチャに現れない**: 出力は module 変数に書かれるため、関数定義だけ見ても「何を返すか」がわからない。ただし caller 側の `use m_foo, only: bar` を見れば出力は特定でき、LLM なら module 宣言部と caller の use only を同時に見て追跡できるため、実質的な問題は小さい
- **module 状態の一括クリーンが面倒**: type なら `deallocate(state)` で一発だが、module 変数は個別に deallocate する subroutine が必要。ただし LLM は module ヘッダの変数宣言を見て漏れなく dealloc subroutine を書けるため、LLM 時代にはこの弱点の実害が縮小している

### 実例: WB.3 リファクタ

`type(wv_storage)` / `type(sxcf_state)` を捨て、module-level singleton 化した結果:
- **行数減**: type 定義、constructor、accessor が不要に
- **OpenACC 制約回避**: structured data region 内の BLOCK 禁止等を自然に回避
- **FILE ↔ MEMORY_3D backend 切替が caller 透過に**: m_wv_storage 内部のフラグ切替だけで、m_sxcf_sc や main_hrcxq は変更不要
- **hgw_combined 統合が容易に**: hrcxq + hsfp0_sc の状態が module 変数として共有できるため、2つのプログラムを1プロセスに統合する際の状態受け渡しが自然
- `comm` のような caller ごとに異なる値は subroutine 引数の optional として直交化

### 例外: m_struc_def のポインタ配列型

Fortran には「allocatable 配列の配列」(ragged array) を直接書く構文がない。
`real(8), allocatable :: a(:)` の配列 `a(n)` は作れない。

`m_struc_def.f90` はこの言語制約を回避するための最小限の wrapper type を定義している:
```fortran
type s_rv1
   real(8),allocatable:: v(:)
end type s_rv1
type s_rv2
   real(8),allocatable:: v(:,:)
end type s_rv2
! ... s_rv4, s_rv5, s_cv3, s_cv4, s_cv5, s_nv2, s_sblock
```

使用例: `type(s_rv2) :: basis(natom)` で、`basis(1)%v` は (10,5)、`basis(2)%v` は (8,3) のように各原子で異なるサイズの配列を持てる (C の `double **a` に相当)。

これは singleton 規約の例外ではなく、**状態を持たない純粋なコンテナ型**。
振る舞い (subroutine) を持たず、module 間のデータフローを隠蔽しない。
singleton が禁じている「状態管理 type」とは明確に区別される。

### 適用ガイドライン
新規モジュールを設計する時:
- まず module 変数 (state) と公開する subroutine を決める
- `type` を導入したくなったら一度立ち止まる — 本当に複数インスタンス必要か?
- caller に `type(...) :: x` を持たせるくらいなら module-level singleton にする
- subroutine 引数は scalar/array of intrinsic types に限定 (state は use で渡す)
- comm のような「caller ごとに異なる値」は subroutine 引数の optional として直交化
- ポインタ配列 (ragged array) が必要な場合のみ m_struc_def のコンテナ型を使用

## ビルド環境 (kt1)

### コンパイラ
- **nvfortran 26.1** (NVIDIA HPC SDK)
- `FC=nvfortran cmake .. -DBUILD_GPU=ON -DBUILD_MP=OFF -DBUILD_MP_GPU=ON`

### 重要な制約
- **`BUILD_MP=ON` は使わない**: rdsigm2.f90 で nvfortran が signal 11 (コンパイラクラッシュ) を起こす。`BUILD_MP=OFF`, `BUILD_MP_GPU=ON` のみ使用
- **nvfortran の間欠的 signal 11**: `-j4` ビルドでランダムにコンパイラがクラッシュする。`make -j4` を数回リトライすれば通る。特定ファイルではなく毎回異なるファイルで起きる
- **cmake の file(GLOB)**: 新規 .f90 ファイル追加時は必ず cmake 再 configure が必要。make だけでは検出されない
- **GPU FP 非決定性**: 同コードで GPU0 と GPU1 の run-to-run で md5 が異なる。回帰比較は `CUDA_VISIBLE_DEVICES` 固定で同条件統一すること

### ディレクトリ構成
```
~/ecaljdeveloper/SRC/                    Fortran ソース
                /SRC/exec/build/         ビルド (GPU=ON, MP=OFF, MP_GPU=ON)
                /SRC/exec/               gwsc, ctrlgenToml.py 等のスクリプト
~/bin/                                   インストール先 (汎用)
~/bin2/                                  GW1500 production 用 (slot_scheduler 対応 run_cmd.py)
```

## GPU 開発の教訓

### WB.4 async 事件 (2026-05-01 → 05-04)

`m_sxcf_sc.f90` に `!$acc async(1)` を追加して GPU kernel を非同期実行する最適化を入れたところ、hsfp0_sc_mp_gpu が生成する sigm が完全に壊れ、全物質の lmf が `zhev_tk4: nev/=nevout` で即死した。

**教訓:**
1. **testecalj PASS は GPU 信頼性の証明にならない** — 2 原子系ではタイミング依存バグが顕在化しない。実物質 (6-8 原子系) でのテスト必須
2. **GPU async は同期ポイントの設計が極めて難しい** — `!$acc wait` が 1 箇所でも不足すると sigm が完全破壊 (中間的劣化ではなく NaN/overflow)
3. **ベンチマーク未実施で production に入ってはいけない** — 性能ゲイン未計測のままリスクを取った
4. **async 再導入時の必須テスト**: GW1500 実物質で sigm の bit-wise 一致テスト、`!$acc wait` の網羅性検証 (全 data region exit 前、全 host_data 前、全 CPU-GPU 境界)

### OpenACC 一般注意
- nvfortran は `!$acc data` (structured) スコープ内の Fortran `BLOCK` を許可しない → `!$acc enter/exit data` (unstructured) で回避
- `attributes(device)` 変数は Fortran スコープで自動解放される → スコープを跨ぐ場合は外部スコープに移動
- `!$acc wait(1)` は goto パス (エラーハンドリング) でも必要

## 入力ファイル: TOML 方式 (2026-05 〜)

### 新方式
Fortran バイナリは以下のみを読む:
- `ctrlG.<sname>.toml` — ctrl + GW driver sections + PB cut-offs の統合 TOML
- `PB.toml` — per-atom product-basis tables (sname-free)

### 生成
```bash
# POSCAR から新規生成
ctrlgenToml.py <sname>        # writes ctrlG.<sname>.toml + PB.toml

# 旧形式 (ctrl + GWinput) からの移行
Legacy2toml.py <sname>        # writes ctrlG.<sname>.toml + PB.toml
```

### 旧 -v オーバーライドの新形式
```
OLD:  lmf si -vnk=8 -vmetal=3
NEW:  lmf si -v[bz.nkabc]=[8,8,8] -v[bz.metal]=3
```

### ファイル命名規約
- `ctrlG.foobar.toml` (foobar = 物質ID, mp-1234 等)
- ディレクトリで識別可能でも**ファイル名に ID を残す** — HTP バッチで取り違え防止、LLM のクロスリファレンス誤認防止

## hgw_combined (in-memory W)

hrcxq (screened interaction W) + hsfp0_sc --job=2 (correlation Sc) を 1 MPI プロセスに統合。
`__WVR/__WVI` ファイル中継 (30-90 GB/物質) を排除し、module-level buffer で in-memory 受け渡し。

gwsc での呼び出し:
```python
# 旧 (3ステップ、ファイル中継)
hsfp0_sc --job=1   # Sx
hrcxq --jobgw=1    # W  →  __WVR/__WVI 書き出し
hsfp0_sc --job=2   # Sc ←  __WVR/__WVI 読み込み

# 新 (1ステップ、in-memory)
hgw_combined --jobgw=1   # Sx + W + Sc (メモリ内)
```

## GW1500 量産インフラ

### スロットスケジューラ
- Unix socket ベース (`/tmp/slot_scheduler.sock`)
- CPU slot ×2 (lmf, np=30 each) + GPU slot ×2
- 6 workers が FIFO 順にスロット取得
- `~/bin2/run_cmd.py` が自動でバイナリ種別を検出し slot request
- `~/bin2/slot_scheduler_daemon.py` が割当管理
- 詳細: `ecalj_auto/README_slot_scheduler.md`

### NaN watchdog
- worker.sh が 5 分ごとに `lsc`, `lsx`, `lqpe`, `llmf` の NaN をチェック
- 8 時間ハードタイムアウト
- 検出時は gwsc + 子プロセスを kill

### 既知の失敗パターン
- NaN in lsc/lqpe: 物質固有の計算困難性 (リトライ不要)
- TIMEOUT 8h: 8 原子系の VRAM 不足等
- bmix minimum: SCF 収束不良
- `/dev/shm/sem.OMPIO*` 残骸: kill 後に残存 → 次の hrcxq がハング。`rm` で即解決
