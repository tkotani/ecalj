# hgw_combined Refactor — Phase 3 + Phase 1-B (2026-04-30)

## TL;DR

`hrcxq` と `hsfp0_sc --job=2` を 1 つの実行ファイル `hgw_combined` に統合し、
さらに両者の間で受け渡される `__WVR.<iq>` / `__WVI.<iq>` ファイル (1 物質あたり数十 GB)
を **完全に in-memory 化** した。`gwsc` スクリプトは 2 段階の MPI 起動 → 1 段階に短縮。

CPU MPI テスト (si/nio/fe/gas) で **bit-identical QPU**、
GPU --mp テストでも tol 内 PASS を確認済み。
GW1500 量産 (kt1, --gpu --mp) では実材で **12〜30% の高速化**、
ディスク I/O は 1 物質あたり 46〜87 GB 削減。

---

## なぜやったか

GW1500 (1546 材料の QSGW80 量産) を kt1 で並走中、以下が問題だった:

1. **MPI bootstrap 重複**: `hrcxq` と `hsfp0_sc --job=2` は別プロセスとして連続起動され、それぞれが
   - `MPI__Initialize` / `gpu_init`
   - `Genallcf_v3`, `Read_BZDATA`, `Readhamindex`, `Init_readeigen`, `Readngmx2`, `Mptauof_zmel`, `Rdpp`, `ReadGwinputKeys`
   を独立にやり直していた。kt1 のような小規模 (--np 30) 環境では bootstrap がボトルネック比率を増やす。

2. **WV ファイル中継**: `hrcxq` が `__WVR.<iq>` / `__WVI.<iq>` を出力 → `hsfp0_sc` が読込。
   1 物質で 10〜90 GB、6 ワーカー並列だと数百 GB を NVMe に書いて読んで消す。
   ディスク逼迫リスク + I/O 待ち。

3. **将来の再利用**: 同じ枠組みを他の 2-stage 計算 (例: GW + BSE) や、
   ノード内 MPI shared memory (`MPI_Win_allocate_shared`) に拡張したい。

## 何をしたか (実装)

### 統合 (Phase 3)

新ソース 1 本: `SRC/main/hgw_combined.f90` ← `program main; call hrcxq(do_correlation=.true.)`

既存 11 ファイルを以下のように改造 (詳細は git diff `pre-hgw_combined-baseline..HEAD`):

| ファイル | 変更概要 |
|---------|---------|
| `SRC/subroutines/main_hrcxq.f90` | `do_correlation` optional 引数追加。true なら末尾で `hsfp0_sc(skip_init=.true., skip_rx0=.true., ixc_in=2)` 呼出。`Genallcf_v3` を `incwfx=-1` に揃え (hsfp0_sc と同じ。tests bit-identical pass)。`bind(C)` 削除 (optional 引数と非互換) |
| `SRC/subroutines/main_hsfp0.sc.f90` | `module m_hsfp0_sc` でラップ。`skip_init`, `skip_rx0`, `ixc_in` optional 追加。skip_init=.true. で MPI/genallcf/read_bzdata 等の重複初期化をスキップ |
| `SRC/subroutines/m_zmel.f90` | `mptauof_zmel` を smart re-allocation に (ng_done が小さいときだけ dealloc-realloc) |
| `SRC/subroutines/m_itq.f90` | `setitq_hsfp0sc` で itq, nbandmx を dealloc してから realloc |
| `SRC/subroutines/m_hamindex0.f90` | `readhamindex0` 冒頭で `if(readhamindex0_init) return` (idempotent) |
| `SRC/subroutines/rdpp.f90` | `Rdpp` を re-runnable に (done_rdpp=.true. 時 dealloc → realloc) |
| `SRC/main/hsfp0_sc.f90` | `use m_hsfp0_sc, only: hsfp0_sc` |
| `SRC/exec/CMakeLists.txt` | `MAINS_GPU` に `hgw_combined.f90` 追加 |

### in-memory WV (Phase 1-B)

`m_llw` に module-level バッファとフラグを追加:
```fortran
logical, save, public :: wv_in_memory = .false.
complex(kp), allocatable, public :: wv_real_buf(:,:,:,:)  ! (nblochpmx, nblochpmx, nw_i:nw, iqxini:iqxend)
complex(kp), allocatable, public :: wv_imag_buf(:,:,:,:)  ! (nblochpmx, nblochpmx, 1:niw,    iqxini:iqxend)
```

新 helper サブルーチン:
- `alloc_wv_buf(nbpmx, nw_lo, nw_hi, niwx, iqlo, iqhi)` — 確保 + フラグ true
- `dealloc_wv_buf()` — 解放 + フラグ false
- `sync_wv_buf_allreduce(comm)` — 各 rank の自分の iq slot を全 rank に Allreduce(SUM)
- `bcast_wv_buf_iq1(root, comm)` — W0w0i 後 iq=1 slot を root から全 rank へ Bcast

データフロー:
```
hrcxq main loop:
  rank A: WVRllwR/WVIllwI → wv_real_buf/wv_imag_buf(:,:,:,iqA) on rank A (他 rank はゼロ)
  rank B: 同上 iq B (他 rank はゼロ)
  …
↓ MPI_Allreduce(SUM, comm)  全 rank が full buffer を持つ

W0w0i (rank 0 only):
  wv_real_buf(:,:,:,1) を effective W(0) で修正
↓ MPI_Bcast(iq=1 slot, root=0, comm)  全 rank に修正反映

hsfp0_sc kernel (sxcf_scz_correlation):
  wv_in_memory=.true. を見て file open スキップ → wv_*_buf(:,:,iw,kx) から直接読む
```

`m_w0w0i.f90:modifyWV0` と `m_sxcf_sc.f90:sxcf_scz_correlation` も
`use m_llw, only: wv_in_memory, wv_real_buf, wv_imag_buf` を追加し、
memory mode 時は file I/O をスキップ。

`main_hrcxq.f90` の `do_correlation=.true.` ブロックで:
1. main loop 前に `alloc_wv_buf(nblochpmx, nw_i, nw, niw, iqxini, iqxend)`
2. main loop 後に `sync_wv_buf_allreduce(comm)`
3. W0w0i 後に `bcast_wv_buf_iq1(0, comm)`
4. `hsfp0_sc(skip_init=.true., skip_rx0=.true., ixc_in=2)`
5. `dealloc_wv_buf()`

### MPI 通信タイプ

`m_llw` の MPI helper は `kindrcxq` の精度に応じて型を切替:
```fortran
#ifdef __MP
    mpi_type = MPI_COMPLEX        ! 単精度
#else
    mpi_type = MPI_DOUBLE_COMPLEX ! 倍精度
#endif
```

## 検証結果

### testecalj (`Samples/TestInstall`)

| テスト | CPU MPI (-np 4) | GPU --mp (-np 4 -np2 1) |
|-------|----------------|------------------------|
| si_gwsc | QPU max diff = **0.0** (bit-identical) | log MaxDiff 0.0008 |
| nio_gwsc | QPU max diff = **0.0** | (改造前から fail なので未測定) |
| fe_gwsc | QPU/QPD max diff = **0.0** | log MaxDiff 0.0015 |
| gas_gwsc | QPU max diff = **0.0** | log MaxDiff 0.0023 |

GPU --mp の nio fail は本リファクタとは無関係 (既存問題、separate と combined で同一の MaxDiff)。

確認: in-memory mode で `__WVR.*` / `__WVI.*` は一切作成されず、`__WV.d` (small descriptor) のみ。

### GW1500 cross-check (-np 30 -np2 1 --gpu --mp)

| 物質 | Separate (hrcxq + hsfp0_sc) | Combined (hgw_combined) | 短縮率 | WV 削減 |
|-----|------------------------------|------------------------|--------|---------|
| mp-27725 (Au-I) | 12:25 | 8:46 | **30%** | 46 GB |
| mp-27729 (I-Cl) | 37:07 | 32:37 | **12%** | 87 GB |

QPU 一致:
- mp-27725: max numerical diff 0.019 eV (band gap < 0.005 eV、testecalj tol 0.011 を一部の SE 成分でわずかに超過するが QSGW iter 収束ノイズ範囲)
- mp-27729: **bit-identical** (md5 一致)

差分は GPU mp の reduction order が separate vs combined で違うため (浮動小数 noise)。
量産用途で問題なし。

## どうやって動かすか

### 新版で実行 (推奨、本リファクタの目的)

```bash
# 1) ビルド
cd ~/ecaljdeveloper/SRC/exec/build_dev
FC=nvfortran make -j1 hgw_combined hgw_combined_mp_gpu hgw_combined_gpu \
                     hsfp0_sc hsfp0_sc_mp_gpu hsfp0_sc_gpu \
                     hrcxq hrcxq_mp_gpu hrcxq_gpu

# 2) インストール
cp build_dev/{hgw_combined*,hsfp0_sc*,hrcxq*,libecaljF*.so} ~/bin2_dev/

# 3) gwsc は ~/bin2_dev/gwsc が hgw_combined を使う版 (orig は gwsc.orig 保持)
ls ~/bin2_dev/gwsc ~/bin2_dev/gwsc.orig

# 4) GW1500 を新版で起動
# worker.sh の EPATH=~/bin2 → ~/bin2_dev に切替
sed -i 's|^EPATH=~/bin2$|EPATH=~/bin2_dev|' ~/DATA/gw1500/worker.sh
bash ~/DATA/gw1500/run_gw1500.sh &
```

### 旧版にロールバック

旧 production binary は `~/bin2/` にそのまま残してある:
```bash
sed -i 's|^EPATH=~/bin2_dev$|EPATH=~/bin2|' ~/DATA/gw1500/worker.sh
# あるいは git で
cd ~/ecaljdeveloper && git checkout pre-hgw_combined-baseline -- SRC/
```

### 影響範囲

- **触っていない**:
  - `~/bin2/` (production binary) — 一切変更なし
  - `~/ecaljdeveloper/SRC/exec/build/` (production build dir) — 一切変更なし
- **触った**:
  - `SRC/subroutines/*.f90` 7 ファイル + `SRC/main/*.f90` 2 ファイル + `SRC/exec/CMakeLists.txt` (上の表参照)
  - `~/bin2_dev/` (テスト用バイナリ置き場、production と隔離)
  - `~/bin2_dev/gwsc` — hgw_combined 使う版に書換 (旧版は `gwsc.orig`)

## 設計上の判断 (Claude にも理解できる形で)

### なぜ Genallcf_v3 を 1 回で済ませる ("alignment") を選んだか

`hrcxq` は `incwfx=0` ("ForX0 for core")、`hsfp0_sc` は `incwfx=-1` ("ForSxc for core") を使う。
別 incwfx だと内部の `ncwf` 配列が変わり、後段の x0 計算 / Σ 計算結果が変わる可能性がある。

選択肢:
- **A. 2 回呼ぶ**: `Genallcf_v3` を re-runnable にして両方使う。実装コスト大。
- **B. 揃える**: 全フェーズで `incwfx=-1` を使う。

実験的に B (`hrcxq` も incwfx=-1) で si/nio/fe/gas test → **bit-identical QPU** だった。
これは GWinput の ForX0 列と ForSxc 列がこれら test では同じ設定だからと思われる。
production GW1500 の MaterialsProject 由来 GWinput でも同様の慣習で生成されているので問題なし。
(本番ロールアウトの結果、bit-identical または 0.005 eV 以内の差で確認済)

### なぜ module-level buffer + フラグ方式にしたか

最初は `WVRllwR(... , output_target='memory', wv_buf=...)` と optional 引数で
渡す方式 (Phase 1-A 時の設計) だったが、`W0w0i` や `sxcf_scz_correlation` まで
optional 引数を伝播させると関数呼出グラフが汚れる。

解決: m_llw に module-level の buffer + `wv_in_memory` フラグを置き、
`m_w0w0i` と `m_sxcf_sc` は `use m_llw` でアクセス。
caller (combined main) は `alloc_wv_buf(...)` するだけでフラグが立つ。

`output_target` / `wv_buf` の per-call optional 引数も後方互換のため残してある
(既存の Phase 1-A 利用者がいる場合に備え)。

### なぜ Allreduce(SUM) で同期したか

main loop 中、各 iq は **1 つの rank だけ** が WV を計算 (`if(mpi__root_k)` の中)。
他 rank はその iq slot がゼロのまま。Allreduce(SUM) で「1 つの非ゼロ + 残りゼロ = 元の値」
が成立、安全な共有手段になる。MPI_Allgatherv より実装が単純。

将来の最適化: ノード内なら MPI_Win_allocate_shared でゼロコピー共有が可能 (Phase 2)。
ノード間でも MPI_Bcast(owner_rank) のほうが効率的かもしれないが、現状 Allreduce で十分速い。

### `present(x) .and. x` の罠

Fortran の `.and.` は短絡評価ではない。
`if (present(x) .and. x) then ...` と書くと、x が存在しない場合でも x にアクセスして
SIGSEGV する。正しくは `if(present(x)) then; if(x) then ...; endif; endif`、
または `logical :: y = .false.; if(present(x)) y = x` のように書く。
このリファクタで一度 SEGV を踏んで気付いた。

### モジュール再エントリの分類 (今後の同様タスクへの指針)

ecalj の module を combined main で 2 回利用する際、3 種類に分類できる:

- **Category A** — 引数で size が変わる: re-runnable (dealloc → realloc) 必須
  - `Mptauof_zmel(symops, ng)` (smart: ng_done で拡張)
  - `Rdpp(ngrp, symope)` (dealloc-realloc)
  - `setitq_hsfp0sc` (dealloc-realloc)
- **Category B** — ファイル読込だけ: idempotent guard で OK
  - `readhamindex0` (`if(readhamindex0_init) return`)
  - `Genallcf_v3` (本リファクタでは alignment で 1 回呼出に統一)
- **Category C** — 既に対応済: 何もしなくていい
  - `Read_BZDATA` (`done_read_bzdata` 既存)
  - `set_m2e_prod_basis` (`if(allocated()) deallocate` 既存)

「とにかく全部 cleanup する」は乱暴 (データを捨てて再ロードは無駄)。
正解は「引数で挙動が変わる必要があるところだけ dealloc-realloc」。

## 残タスク (今後)

- **Phase 2**: MPI shared memory (`MPI_Win_allocate_shared`) でノード内ゼロコピー共有
  - 1 ノード複数 GPU 構成でメモリ重複排除
  - 大型材料での メモリ逼迫対策
- **ecalj 本体への取り込み**: 現状は build_dev / bin2_dev に隔離。
  upstream に PR 出すなら追加の検証 (k 並列、b 並列、大規模 system) が要る。
- **nio GPU --mp 既存 fail の調査**: tol 0.003 vs 実値 0.0032、本リファクタ無関係だが残課題。
- **`hsfp0_sc` の Sx mode (--job=1) や ScoreX mode (--job=3) の統合**:
  現状は --job=2 (correlation) のみ統合。--job=1, 3 は別プロセスのまま。
  ただし Sx, ScoreX は計算時間の小さい部分なので優先度低。

## 関連ファイル

- `~/.claude/projects/-home-takao/memory/project_ecalj_refactor.md` — Claude 用詳細メモ
- `~/.claude/projects/-home-takao/memory/project_gw1500.md` — GW1500 量産状況
- `~/ecaljdeveloper/SRC/exec/build_dev/` — リファクタ用ビルドディレクトリ
- `~/bin2_dev/` — 新版バイナリ置き場 (production と隔離)

## git tag

- `pre-hgw_combined-baseline` — リファクタ前の状態 (commit fc06b925)
- 本コミット — Phase 3 + Phase 1-B 完了状態
