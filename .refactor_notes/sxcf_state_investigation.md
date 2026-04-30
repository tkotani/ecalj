# sxcf_state リファクタリング 調査レポート (Gate 1)

調査日: 2026-04-30 (m_itq refactor commit `5cebb798` 後の状態)

ゴール: m_sxcf_sc.f90 の module-level state を整理して streaming Phase 1-C の再挑戦を可能にする (= [2])。
本ノートは「どこを安全に/危険に触れるか」の事前マッピング。

---

## 1. OpenACC directive のかかっている変数

`grep "!$acc" m_sxcf_sc.f90` の結果:

| 変数 | OpenACC 操作 | 場所 | 移設リスク |
|---|---|---|---|
| `zsecall` | `enter/exit data create`, `kernels`, `host_data use_device` | sxcf_scz_exchange (line ~167-247), sxcf_correlation_init (~312), step_kx (~628), finalize | **高** |
| `wvi_upper` | `enter data copyin`, `present`, `exit data delete` | step_kx (SetWVblock, iwimag, ReleaseWV) | **高** |
| `wvr_upper` | 同 | 同 (iwreal) | **高** |
| `idx_i, idx_j` | 同 | step_kx | **高** |
| `zmel` | `host_data use_device`, `present` (m_zmel が所有) | step_kx, exchange | n/a (m_zmel) |
| `wgtim, wgtiw, nttp, itw, itpw` | `data copyin` | step_kx (内部 block 内ローカル) | n/a (block local) |
| `vcoud, wklm(1), wk(1), wtff` | `data copyin` (exchange) | sxcf_scz_exchange | n/a (m_readVcoud / m_read_bzdata) |

**OpenACC を触っていない変数 = 移設「安全」**:
- `ekc, eq, omega` (real(8), allocatable)
- `t_sw_zmel, t_sw_xc, t_sw_cr, t_sw_ci, t_sw_setwv` (type(stopwatch))
- `emptyrun, keepwv` (logical scalar)
- `kx, irot, ip, isp, ntqxx, nt0p, nt0m, ifrcw, ifrcwi, wkkr` (loop iterators / per-iteration scalars)
- `ndiv, nstatei, nstatee` (allocatable, currently set/used in init only)

## 2. sxcf_scz_exchange と sxcf_scz_correlation の共有 module 変数

各 subroutine 内での参照回数 (declaration や comment は除く実使用):

| 変数 | exchange | correlation (init+step_kx+finalize) | 共有? |
|---|---|---|---|
| kx | 12 | 21 | ✓ |
| irot | 10 | 10 | ✓ (loop) |
| ip | 12 | 13 | ✓ (loop) |
| isp | 16 | 16 | ✓ (loop) |
| ntqxx | 5 | 13 | ✓ |
| nt0p, nt0m | 0 | 4, 3 | correlation only |
| ifrcw, ifrcwi | 0 | 6, 4 | correlation only |
| wkkr | 2 | 3 | ✓ |
| ekc | 4 | 9 | ✓ |
| eq | 3 | 3 | ✓ |
| omega | 0 | 12 | correlation only |
| emptyrun | 1 | 1 | ✓ |
| keepwv | 0 | 4 | correlation only |
| wvr_upper, wvi_upper | 0 | 10, 10 | correlation only |
| t_sw_zmel, t_sw_xc | 5, 5 | 5, 5 | ✓ |
| t_sw_cr, t_sw_ci, t_sw_setwv | 0 | 5, 5, 11 | correlation only |
| zsecall | 8 | 9 | ✓ |
| ndiv, nstatei, nstatee | 0 | (set in init's PreIcountBlock, freed in same block) | local-ish |

**結論**: 共有変数を derived type に移すには exchange と correlation を**同時に**書き換える必要。correlation 専用変数は片方だけ書き換えても安全。

## 3. 触りやすさのクラス分け

### Class A: 移設しやすい (correlation 専用 + OpenACC なし)
- `omega` (real, allocatable)
- `keepwv, debug` (logical)
- `nt0p, nt0m, ifrcw, ifrcwi` (integer scalar)
- `t_sw_cr, t_sw_ci, t_sw_setwv` (stopwatch)

これらは sxcf_correlation_init / step_kx / finalize の内部だけで完結する。derived type に
入れたり、step_kx ローカルに降ろしたりが**他の関数を一切壊さない**。

### Class B: 移設しやすいが共有 (両 subroutine で使う + OpenACC なし)
- `ekc, eq` (real, allocatable)
- `ntqxx, wkkr` (scalar)
- `t_sw_zmel, t_sw_xc` (stopwatch)
- `emptyrun` (logical)
- `kx, irot, ip, isp` (loop iterators)

derived type に入れるなら exchange と correlation の両方を更新する必要。
Loop iterators は両方独立でローカル化可能 (両方 do kx=1,nqibz の自前ループがあるので、各々
ローカル変数 `kx_loc` にして使えば両者独立で動く)。

### Class C: 移設リスク高い (OpenACC directive 多数)
- `zsecall` (4D complex, output)
- `wvi_upper, wvr_upper, idx_i, idx_j` (correlation 専用、step_kx 内で alloc/dealloc)

`!$acc enter data create(state%ws%zsecall)` のような derived type member への directive は
nvfortran で動く保証がない。**この調査では検証していない**。
保守的には:
- `zsecall` は module-level のまま
- `wvi_upper, wvr_upper, idx_i, idx_j` は step_kx ローカル化に留める (block-local の openacc enter/exit data)

## 4. block 構文の host scope 確認

m_sxcf_sc 内で `kx, irot, ip, isp` を subroutine ローカルに移しても、`SetWVblock`,
`get_correlation_block`, `ReleaseWV`, `CorrelationSelfEnergyImagAxis`,
`CorrelationSelfEnergyRealAxis` などの inner `block` 構文は、Fortran 仕様上 **host scope
の変数は見える** (block は「サブプログラム」ではなく「式の集まり」)。

ただし inner block が**同名の変数を declare**するとシャドーする。grep で確認:

- inner `iqini, iqend` (SetWVblock): 内部 declare → そのスコープで上書き OK
- inner `iw, i, j, tri_idx, ittp` (各種 block): 内部 declare → 局所
- `ns1, ns2, kr, ip, isp, irot, kx`: inner block で declare されていない → host scope を見る

**結論**: kx/irot/ip/isp を subroutine ローカル化しても安全。

## 5. snapshot baseline 比較スクリプト

`~/ecaljdeveloper/.refactor_notes/snapshot_baseline.sh capture <label>` で、
testecalj を走らせた直後の各 *_gwsc_work 内に対し:

1. `SEBK/SECU, SEC2U, SEXU, SEX2U` の md5
2. `lcombined` の主要診断行 (ntq=, nbandmx, ncount, ω 分割)
3. `SEBK/SECU` の最初 100 行 (浮動小数 D24.16 形式)
4. `lcombined` の末尾 50 行

を記録する。`diff <a> <b>` で 2 つの snapshot を比較。

baseline (m_itq refactor 後) は `pre_step2` ラベルで取得済み:
- `~/ecaljdeveloper/.refactor_notes/snapshots/pre_step2/{si,nio,fe,gas}_gwsc_work.txt`

## 6. 推奨ステップ (根拠つき)

> 各ステップ後 `snapshot diff pre_step2 step2.X` で内部状態の bit-equal 確認。

### Step 2.1 — Loop iterator ローカル化 (Class B 部分)
変数: `kx, irot, ip, isp` を `sxcf_scz_correlation` および `sxcf_scz_exchange` のローカル変数に降ろす。module 宣言から削除。
- 効果: 「いま module の kx は何の値だ?」混乱が消える
- 規模: 各 subroutine で `do kx=1,nqibz` を `do kx_loc=1,nqibz; kx = kx_loc` のように外側で代入してから使う、もしくは Body の `kx` を `kx_loc` に置換 (sed で大量置換)
- リスク: low. block 構文も host scope で見えるので壊れない
- テスト: snapshot diff bit-equal

### Step 2.2 — stopwatch の type 化 (Class A + B)
5 個の `t_sw_*` を `type :: sxcf_stopwatches` にまとめる。最初は **module-level** に `type(sxcf_stopwatches) :: sw` として置く (まだ derived type を引数で渡すアプローチには行かない)。`t_sw_zmel` を `sw%zmel` にリネーム、等。

- 効果: derived type を試運転、5 個の global がまとまる
- リスク: low. stopwatch は表示時間以外に副作用ないし、両 subroutine が同じ `sw%zmel` を見れば現状と等価
- テスト: snapshot diff bit-equal

### Step 2.3 — correlation 専用変数の workspace type 化 (Class A)
`omega, keepwv, nt0p, nt0m, ifrcw, ifrcwi, t_sw_cr, t_sw_ci, t_sw_setwv` を
`type :: sxcf_correlation_workspace` にまとめる。**module-level に保持**するが struct で
グルーピングするだけ。step_kx 内部の参照を `ws%omega` 形式に。

- 効果: correlation だけ綺麗になる (exchange に影響なし)
- リスク: low. exchange は触らない
- テスト: snapshot diff bit-equal

### Step 2.4 — Class B 共有変数の grouping (任意)
`ekc, eq, ntqxx, wkkr, emptyrun` も derived type に入れる。**両 subroutine を更新する必要**あり。
exchange も `ekc → ws%ekc` に書き換え。

- 効果: 共有 state がスコープ可能になる
- リスク: medium. 両方更新する一括変更
- テスト: snapshot diff bit-equal (両 SEX/SEC ファイル)

### Step 2.5 — derived type を引数化 (state machine 化、本命)
ここまでは module-level の type instance だった。これを `subroutine sxcf_correlation_init(state, ...)` のように引数渡しに切り替える。
exchange 用も `subroutine sxcf_exchange(state, ...)`。streaming で複数 state を持てる構造になる (実利用の予定はないがテスタブル)。

- 効果: state machine 化の本命。streaming Phase 1-C 再挑戦の地盤完成
- リスク: high. 全ての caller を更新。intent semantics 注意
- テスト: snapshot diff bit-equal、Phase 1-C 再挑戦

### Step 2.6 — `zsecall` の所有検討 (任意・後で)
Class C なので最後。`zsecall` を state 内に入れるかは別判断 (OpenACC 互換性次第)。

## 7. 2.4 / 2.5 を「やらない」基準

* 2.4 で snapshot diff が出たら、変数の正確な意味論が共有変数間で違う可能性 → 深堀り必要
* 2.5 で intent / 寿命周りの罠が見えたら、type 引数化を諦めて 2.4 までで止める

## 8. 不確定事項 (次の調査になったら確認したい)

- nvfortran が `!$acc enter data create(state%ws%array)` を許すか (Class C 移行に必要)
- `MPI__reduceSum` が derived type の component 引数でも動くか (ABI 上は普通の配列なので OK と推測)
- `protected, target` 属性付きの allocatable component が derived type で許されるか (Fortran 標準上 OK)
