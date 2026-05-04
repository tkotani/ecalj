# ecalj Claude Development Guide

ecalj を Claude (LLM) で開発する際のガイドライン。
コーディング規約、ビルドの罠、GPU 開発の教訓、アーキテクチャ方針をまとめる。

## コーディング規約: module = singleton (kotani 設計)

ecalj の Fortran コードは **module-as-singleton-class** パターンで書く。

### 原則
- **1 module = 1 責務 = 1 状態セット**: 複数インスタンス不要なら `type` を被せない
- **状態は module 変数で公開**: `protected, public` で読み取り公開、書き換えは module 内 subroutine 経由のみ
- **import は明示**: 必ず `use m_foo, only: bar, baz` で必要な symbol だけ取る
- **subroutine 規約**: 入力=引数、出力=module 変数で返す
- **`type` (derived type) は避ける**: Fortran の type には以下の罠がある
  - `intent(out)` で allocatable components が暗黙に deallocate される
  - `present(state%alloc_member)` が nvfortran + OpenACC で silent fail
  - 合成 type の copy semantics が不明瞭
  - `use only` でデータフローが見えない

### 実例
WB.3 リファクタで `type(wv_storage)` / `type(sxcf_state)` を捨て、module-level singleton 化:
- 行数減、OpenACC 制約回避
- FILE ↔ MEMORY_3D backend 切替が caller 透過に
- `comm` のような caller ごとに異なる値は subroutine 引数の optional で直交化

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
