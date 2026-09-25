# MLO-QSGW の理論確認（2026-09-25）

**図**（クリックで VSCode が開きます）:
[mlo_rows.png](Samples/kBT/LiTi2O4/mlo_rows.png) 従来 vs MLO の縦比較 ・
[frozen_vs_not.png](Samples/kBT/LiTi2O4/frozen_vs_not.png) χ̃ 凍結の有無 ・
[sigr_decay.png](Samples/kBT/LiTi2O4/sigr_decay.png) Σ(R) の減衰 ・
[three_routes.png](Samples/kBT/LiTi2O4/three_routes.png) 3 経路 × iter1/2 ・
[mlo_iter_grid.png](Samples/kBT/LiTi2O4/mlo_iter_grid.png) 旧版の 10 反復
（詳細は [§7 図](#7-図)）

## 1. PMT 基底のセットアップ

そのとおりです。MT 内の**球対称**ポテンシャルから $\phi_l(r,E_\nu)$, $\dot\phi_l$ を作り、
それで smooth Hankel 包絡（MTO）を augment し、APW を足す。
基底は球対称 DFT ポテンシャルだけで決まります。

ただし 2 点:

- **密度が変われば基底も変わる**（弱いが変わる）。
- `pwmode = 11` では APW の集合が $q$ に依存する（$|q+G|$ カット）。

## 2. Σ 行列の保持形 — MLO と従来で形が違う

### MLO 経路

$$\Sigma^{\rm MLO}_{\alpha\beta} = c^\dagger \Sigma^\psi c,
\qquad c = \langle\psi|\tilde\chi\rangle$$

なので、**文字どおり $\langle\tilde\chi_\alpha|\hat\Sigma|\tilde\chi_\beta\rangle$**（行列要素）です。
戻すときは $A = S^{\rm PMT} z^{\rm MLO} = \langle\chi|\tilde\chi\rangle$ を使って

$$\mathrm{senex} = A\,O^{-1}\,\Sigma^{\rm MLO}\,O^{-1}A^\dagger
= \langle\chi|\hat P\hat\Sigma\hat P|\chi\rangle$$

（設計書 式 (8)(14)）。整合しています。

### 従来 `sigm`

MTO ブロックのみ（`ndimsig = nlmto`、APW は入らない）で、しかも**両側に $O^{-1}$ を畳み込んだ形**。
だから今朝オンサイトが 73.8 eV もありました（MLO 側は 187 meV）。
**同じ「⟨基底|Σ|基底⟩」ではありません。**

## 3. ⟨ψ|Σ|ψ⟩ から基底表現へ

そのとおりです。ただし窓 $\bar\theta$ の外の状態は落とされ、残りは `eseavrmean` の定数で埋めるので、
**変換が厳密なのは窓の中だけ**です。

---

# 4. いつ MLO を作るのか — 1 反復のタイムライン

## 4.0 用語

| | 中身 | 決まる時期 |
|---|---|---|
| **MLO 索引** | `ix`（どの MTO チャネルを種にするか）、`ndimMTO`、窓のパラメータ `eferm`/`ecbot`/`nskip` | **連鎖の最初に 1 度だけ**（`HamRsMLO`）。$\Sigma^{\rm MLO}$ はこのチャネル上の行列なので、`ndimMTO` を途中で変えられない |
| **MLO 係数** $z^{\rm MLO}(q)$ | $\tilde\chi_\alpha(q)=\sum_m\chi^{\rm PMT}_m(q)z_{m\alpha}$ の展開係数 | **毎反復作り直す**（下記 ①） |

以降 $z_n$ = 反復 $n$ で作った係数、$\Sigma_n$ = 反復 $n$ で作った $\Sigma^{\rm MLO}$。

## 4.1 連鎖の最初に 1 度だけ

| 段 | 実行 | 出力 |
|---|---|---|
| 0a | `lmf --jobgw=0` | `HAMindex0`（構造のみ） |
| 0b | `qg4gw --job=1` | `QGpsi`（offset-Γ 込みの q リスト） |
| 0c | `lmf --writeham --mlo` | `HamiltonianPMTInfo` |
| 0d | `mlo --mlo` | **`HamRsMLO`** — 末尾レコードに MLO 索引（`ix`, `ndimMTO`, `mlomethod`, `nskip`, `eferm`, `ecbot`） |

**ここで決まるのは索引だけ**で、$z^{\rm MLO}$ は作りません。

## 4.2 反復 $n$ の中の順序 ← ここが肝心

`gwsc 1` の 1 回がこの順に走ります。

| # | 実行 | $\tilde\chi$ に関して何が起きるか |
|---|---|---|
| ① | **`lmf --jobgw=1`**（[sugw.f90](SRC/subroutines/sugw.f90)） | |
| | ├ 各 $q$ で $H^{\rm LDA}(q)$, $S(q)$ を作る | |
| | ├ **`getsenex`**: `SigRsMLO`（$=\Sigma_{n-1}$）を読む | **`ZmloSig` の $z_{n-1}$ を使う**。$\Sigma_{n-1}$ はこの基底で書かれたので、これでしか正しく戻せない |
| | ├ 対角化 → $\psi$（GW に渡る波動関数） | |
| | └ **段 a'** | **ここで $z_n$ を作る。** その反復の $H(q)$ から `Hreduction` で構成し、`ZmloNew.<rank>` に書き出す。$c^{\rm MLO}=\psi^\dagger S z_n$ を `__cmlo.data` へ。**$H$ は $\Sigma$ 入り**（下記）|
| ②a | `heftet --job=1` | $E_F$（四面体法）。MLO には無関係 |
| ②b | `hbasfp0 --job=3` → `hvccfp0 --job=3` → `hsfp0_sc --job=3` | コアの交換項 $\Sigma_x^{\rm core}$。**①の $\psi$ を使う**。MLO には無関係 |
| ②c | `hbasfp0 --job=0` → `hvccfp0 --job=0` | 価電子の積基底と Coulomb 行列 $v$ |
| ②d | **`hgw --jobgw=1`** | ①の $\psi$ から $\chi_0 \to W \to \Sigma_x, \Sigma_c$。出力は**バンド添字の $\Sigma^\psi_{ij}(q)$**（`SEXU`, `SEXcoreU`, `SEC*`）。**まだ MLO は関係しない** |
| ③ | **`hqpe_sc`** | ここで初めて MLO に落ちる。`__cmlo.data` の $c_n$ を読み、$\Sigma^{\rm MLO}_n = c_n^\dagger\Sigma^\psi c_n$ → `__SigmMLO.q`。**$z_n$ の基底**。同時に従来の `sigm`（MTO）も並行して書く |
| ④ | **`mlo --mlofreeze --mlo`** | FFT して `SigRsMLO` を書く。**$z_n$ の基底** |
| ⑤ | **昇格**（`gwsc` が rename） | `ZmloNew.*` → `ZmloSig.*`。**この瞬間に「現行の $\tilde\chi$」が $z_{n-1}$ から $z_n$ に切り替わる** |
| ⑥ | **`lmf`**（Σ 入り SCF） | `getsenex` が `SigRsMLO`（$=\Sigma_n$）を **$z_n$** で読む |
| ⑦ | **`job_band ... --mlo`**（掘ったディレクトリの中で） | 同じく $z_n$。ただし Σ メッシュ外の $k$ は `ZmloSig` に無いので、その場の $H$ から作る（読み戻す相手が無いため） |

### ⑦の注意 — 2 つとも今日ハマった

**(1) `--mlo` を必ず渡す。** `job_band` は `parse_known_args()` で残りの引数をそのまま `lmf` に
渡します（[job_band:91,93](SRC/exec/job_band#L91-L93)）。`--mlo` が無いと `c0_mlo` が偽になり、`sigmlo_init` は
冒頭の `if(.not.c0_mlo) return` で即座に戻ります。**MLO 経路は完全に無効**になり、`getsenex` は
従来の `sigm` を使う。`hqpe_sc` は `sigm` と `__SigmMLO.q` を並行して書くので `sigm` は常に
存在し、**黙って従来内挿のバンドが描かれます**。今日の図はすべてこれでした。

**(2) 連鎖ディレクトリの中で走らせない。** `job_band` の第 1 段は `lmf --quit=band` で、
これは **`efermi.lmf` を書き換えます**（[m_bndfp.f90:274](SRC/subroutines/m_bndfp.f90#L274)）。MLO の窓は `efermi.lmf` から
`ecbot` オフセットを読むので、連鎖ディレクトリで描くと窓の基準が動きます。
`snap/iter<N>/` を掘って状態一式を置き、その中で走らせる。**そのままスナップショットになります。**

**答え**: MLO（係数 $z$）は **反復ごとに 1 回、①の段 a' で、その反復のハミルトニアンから**作られます。
そして**⑤で初めて「現行」になります**。

### どの $H$ から作るか — $H^{\rm LDA}$ は初段だけ

$\tilde\chi$ は、これから $\Sigma$ を表現する相手のハミルトニアンの多様体を張っていなければ意味が
ありません。反復 2 以降の $H$ は既に $\Sigma$ を含んでいるので、**$z_n$ も $\Sigma$ 入りの $H$ から作ります**。

| 反復 | a' が使う $H$ |
|---|---|
| 1 | $H^{\rm LDA}$ — まだ $\Sigma$ が無いので（`sigmamode` が偽） |
| 2 以降 | $H = H^{\rm LDA} + \mathrm{senex}(\Sigma_{n-1})$ |

コード上は `zhev_tk4` が `hamm` を破壊するため、senex を足した直後の $H$ を `hamm_qsgw` に
退避して a' で使います（[sugw.f90](SRC/subroutines/sugw.f90)）。2 スロット化しているので、`getsenex` 側は `ZmloSig` を
読むだけであり、a' がどの $H$ を使っても読み戻しは壊れません。

**これは標準の MLO 手法と同じです。** [m_bandcal.f90](SRC/subroutines/m_bandcal.f90) では

| 行 | 処理 |
|---|---|
| [207](SRC/subroutines/m_bandcal.f90#L207) | `hambl` → $H^{\rm LDA}$ |
| [210–211](SRC/subroutines/m_bandcal.f90#L210-L211) | `sigmamode` なら `getsenex` → **`hamm = hamm + senex`** |
| [221〜](SRC/subroutines/m_bandcal.f90#L221) | `if(writeham)` → **`__HamiltonianPMT` に書く** |

の順で、**`__HamiltonianPMT` には senex を足した後の QSGW ハミルトニアンが入ります**。
`mlo` はこれを読んで `HamRsMLO` を作るので、Si の QSGW バンドなどで使っている標準の MLO は
**元から QSGW の $H$ で作られています**。段 a' が `hamm_lda` を使っていたのが標準から外れていた
方であり、今回の変更はそこを揃えたことになります。

### なぜ①の `getsenex` が前反復の $z_{n-1}$ を使うのか

**前反復で作った MLO から、$c$ も $\Sigma$ も導かれているから**です。反復 $n-1$ の中で:

$$z_{n-1}\ \xrightarrow{\ \text{段 a'}\ }\ c_{n-1}=\psi^\dagger S\,z_{n-1}\ \xrightarrow{\ \texttt{hqpe\_sc}\ }\ \Sigma_{n-1}=c_{n-1}^\dagger\,\Sigma^\psi\,c_{n-1}$$

と一本の鎖になっています。$c_{n-1}$ は $z_{n-1}$ から作られ、$\Sigma_{n-1}$ はその $c_{n-1}$ から作られる。
つまり **`SigRsMLO` の数値は $z_{n-1}$ の刻印を持っています** — 中身は
$\langle\tilde\chi^{(n-1)}_\alpha|\hat\Sigma|\tilde\chi^{(n-1)}_\beta\rangle$ であって、
$\tilde\chi^{(n-1)}$ が何者かを知らなければ意味を持ちません。

だから反復 $n$ の①でこれを PMT 基底へ戻すときは、**必ず同じ $z_{n-1}$** で

$$\mathrm{senex}=A_{n-1}\,O_{n-1}^{-1}\,\Sigma_{n-1}\,O_{n-1}^{-1}A_{n-1}^\dagger,
\qquad A_{n-1}=S\,z_{n-1}$$

としなければなりません。ここに $z_n$ を使うと、$\Sigma_{n-1}$ の数値を**別の関数の間の行列要素**として
読み替えることになり、演算子そのものが変わってしまいます。**これが今日の破綻の正体でした。**

一方、段 a' がこの同じ①の中で $z_n$ を作るのは、**次の** $\Sigma_n$ をその基底で書くためです。
だから①の時点で $z_{n-1}$（読む用）と $z_n$（書く用）が同時に要る。

```
    反復 n-1                    反復 n                      反復 n+1
       |                           |                            |
   Σ_{n-1} を z_{n-1} で書く  →  ① getsenex: z_{n-1} で読む      |
                                 ① a'      : z_n を作る          |
                                 ③④      : Σ_n を z_n で書く    |
                                 ⑤ 昇格   : 現行 = z_n          |
                                 ⑥ lmf    : z_n で読む      →  ① getsenex: z_n で読む
                                                                ① a' : z_{n+1} を作る
```

$\tilde\chi$ は**ハミルトニアンに追随して毎反復更新される**が、
$\Sigma$ は**書かれた基底以外で読み戻されることが無い**。両立しています。

## 4.3 なぜ 2 スロット要るのか

①の中で $z_{n-1}$（読む用）と $z_n$（書く用）が**同時に生きている**必要があるからです。
1 つしか持てないと、どちらかが必ず食い違います。

| ファイル | 中身 |
|---|---|
| `ZmloSig.<procid>` | 現行 `SigRsMLO` が書かれた $\tilde\chi$。[m_sigmlo.f90 `zmlo_sig_load`](SRC/subroutines/m_sigmlo.f90) が読む。追記しない |
| `ZmloNew.<procid>` | 段 a' がこの反復の $H$ から作った $\tilde\chi$（[m_sigmlo.f90 `zmlo_new_append`](SRC/subroutines/m_sigmlo.f90)）。⑤で [gwsc](SRC/exec/gwsc) が rename |

ランクごとのファイルなのは、段 a' が MPI で $q$ を分担するのと、
`pwmode=11` で `ndimh` が $q$ ごとに変わり固定レコード長が使えないためです。

## 4.4 コード上どこで MLO を作るか

**作る実体は 1 箇所**、[`Hreduction`](SRC/subroutines/m_hreduction.f90#L45)（`m_hreduction.f90`）です。
`zMLO` を返します。これを 3 種類の場所から呼んでいます。

| 呼ぶ場所 | いつ | 何を決めるか |
|---|---|---|
| [m_HamPMT.f90](SRC/subroutines/m_HamPMT.f90) | 段 0d（`mlo --mlo`）、および MLO バンド | **索引**: `ix`（どの MTO チャネルを種にするか）、`ndimMTO`、窓 `nskip`/`eferm`/`ecbot`。`HamRsMLO` の末尾レコードに保存 |
| [sugw.f90](SRC/subroutines/sugw.f90) 段 a'（`NewChiForThisIteration`）| 反復ごと、`lmf --jobgw=1` の中、各 $q$ で | **書く側の $z_n(q)$**。反復 2 以降は `hamm_qsgw`（Σ 入り $H$）、初段は `hamm_lda` から。`ZmloNew.<procid>` へ |
| [m_sigmlo.f90](SRC/subroutines/m_sigmlo.f90) `getsenex`（`ChiTildeCache`）| `ZmloSig` に無い $k$ のとき | **読む側の $z(k)$**。その場の $H$ から作り、その `lmf` 実行中だけメモリに保持 |

つまり:

- **索引（どのチャネルか・窓）は連鎖の最初に 1 回**だけ決まり、以後凍結（`--mlofreeze`）。
  $\Sigma^{\rm MLO}$ はそのチャネル上の行列なので `ndimMTO` を途中で変えられない。
- **係数 $z^{\rm MLO}(k)$ は各 $k$ で毎回作る**。$\tilde\chi$ は「保存された表」ではなく
  **$H(k)$ に処方を当てる関数**である。だから任意の $k$（バンドプロット）でも作れる。
- `ZmloSig`/`ZmloNew` は「その処方を**いつの $H$** で当てたか」を固定するための置き場であって、
  $\tilde\chi$ そのものの定義ではない。

**$\tilde\chi$ を規定する材料**（`Hreduction` の入力）は

1. $H(k)$, $S^{\rm PMT}(k)$ — その時点のハミルトニアンと重なり
2. `ix` — 種にする MTO チャネル（凍結）
3. 窓 $\bar\theta$ のパラメータ `eferm`, `ecbot`, `nskip`, `mlo_w`, `mlo_delta`, `mlo_wfrz`（凍結）

の 3 つだけです。1 が反復ごとに動き、2・3 は固定、という構図。

## 4.5 MLO 法の要点 — 入力・作り方・内挿の成立条件

### user のまとめ（そのとおり）

1. 入力として**メッシュ k 点での $H^{\rm PMT}(k)$**（と $S^{\rm PMT}(k)$）が要る。
2. それがあれば $H^{\rm MTO}(k)$ は部分ブロックとして取り出せる。いま考えている MTO チャネル
   （索引 `ix`）だけを使って、メッシュ k 点で $|{\rm MLO}\rangle$ を $|{\rm PMT}\rangle$ から作る行列
   $z^{\rm MLO}(k)$ が得られる。
3. $\langle {\rm PMT}|X|{\rm PMT}\rangle$ があれば $\langle {\rm MLO}|X|{\rm MLO}\rangle$ に展開でき
   （部分空間にはなる）、これなら FFT で任意の $k$ に内挿できる。$H$ も $S$ も $\Sigma$ も。

コードと突き合わせると 1・2 は完全に一致する。3 も原理としては正しいが、**成立条件が一つある**。

### 位相（ゲージ）の心配は不要

[m_hreduction.f90:322](SRC/subroutines/m_hreduction.f90#L322):

```fortran
cmlo_loc = matmul(Amat, matmul(evecmto^dagger, ovlmx(ix,ix)))
```

コメントどおり $|F^{\rm MLO}_k\rangle = \hat P\,|F^{\rm MTO}_k\rangle$ であり、
$\Psi^{\rm MTO}_j$ も $\Psi^{\rm PMT}_i$ も **$|\Psi\rangle\langle\Psi|$ の形でしか現れない**ので、
対角化の任意位相は相殺する。種は**裸の MTO 基底関数** $\chi^{\rm MTO}_k$（k 非依存の実空間関数の
ブロッホ和）。したがって $\tilde\chi$ にゲージ依存性は無い。

### 成立条件: $\hat P(k)$ の k についての滑らかさ

$$\hat P(k)=\sum_{ij}|\Psi^{\rm PMT}_i(k)\rangle\,\bar\theta_{ij}(k)\,
\langle\Psi^{\rm PMT}_i(k)|\Psi^{\rm MTO}_j(k)\rangle\langle\Psi^{\rm MTO}_j(k)|$$

$\tilde\chi(k)=\hat P(k)\,\chi^{\rm MTO}(k)$ は、「k 非依存の関数のブロッホ和」を
**k ごとに違う射影で修正したもの**である。q 周期性は保たれるので FFT 自体は正当だが、

> **実空間で短距離になるかどうかは、$\hat P(k)$ が k についてどれだけ滑らかかで決まる。**

では $\hat P$ の k 依存性はどこから来るか。[m_hreduction.f90:308](SRC/subroutines/m_hreduction.f90#L308):

```fortran
Amat(i,j) = fac(i,j) * max( fermidist((evl(i)-efrz)/ewfrz), fermidist((evl(i)-ecut)/ewuse) )
```

**窓ではない。** 数字で確認した（user 指摘）:

- **第 1 項は恒等的に 0**。`mlomethod=4` では [`efrz = -1d99`](SRC/subroutines/m_hreduction.f90#L246)
  （コメント「hard-freeze edge; active only for mlomethod=3」）。効くのは**第 2 項だけ**。
- 第 2 項の位置は $\varepsilon^{\rm cut}_j=\max(\varepsilon_{\rm cbot}+\Delta,\ \varepsilon^{\rm MTO}_j)$。
  **MTO のみの固有値は PMT のそれよりそれなりに上にある**ので、通常は $\varepsilon^{\rm MTO}_j$ が選ばれる。
  するとシグモイドの引数は $(\varepsilon^{\rm PMT}_i-\varepsilon^{\rm MTO}_j)/w$ という**バンド差**になり、
  両者が一緒に分散するぶん k について滑らかになる。
- 幅も狭くない。LiTi₂O₄ の設定は **`mlo_w = 2.0` eV、`mlo_delta = 2.0` eV**（`ctrlg` の
  `[gw]` 節）で、t2g の帯幅 0.4 eV の **5 倍**。

**したがって $\hat P(k)$ の k 依存性は窓由来ではなく、
固有ベクトル $|\Psi^{\rm PMT}_i(k)\rangle\langle\Psi^{\rm PMT}_i(k)|$ 自体が k で回ることから来る。**
これは band 多様体の性質であって、窓のパラメータで緩和できるものではない。

（2026-09-25 の初稿では「窓が狭いのではないか、広げてみては」と書いたが、
`mlo_w = 2.0` eV を確認して撤回した。）

## 4.6 自己無撞着にどう持っていくか — 採る方式

### 出発点: $\langle\chi^{\rm PMT}|\tilde\chi\rangle$ は内挿できない

$z^{\rm MLO}$ の入った行列は **PMT 基底の行**を持つ。設計書 §3.4 のとおり、PMT のうち
q 周期的なのは **MTO 部分だけ**で、**APW 部分は q 周期的でない**
（`pwmode=11` は $|q+G|$ カットなので G の集合自体が q で変わる）。
行が k ごとに別の関数を指すので、**この行列は実空間に落とせず、内挿できない**。

同じ理由で $H^{\rm PMT}(R)$ を保存して任意の $k$ で $H^{\rm PMT}(k)$ を再構成する道も無い。
**$H^{\rm PMT}(k)$ は、その密度で実際にハミルトニアンを組む以外に入手できない。**

→ したがって **$\tilde\chi$ は毎回作るしかない。そして作るにはハミルトニアンが要る。**

### 採る方式

> **反復 $n$ で $\Sigma_n$ を作るとき、その場の $H$（＝反復 $n$ 開始時点のハミルトニアン）から
> $z^{\rm MLO}(k)$ を「後で読むことになる全ての $k$」で作り、保存する。
> 以後 $\Sigma_n$ を読むときは必ずその $z$ を使う。**

$\Sigma^{\rm MLO}(R)$ 自体は MLO 基底上の行列なので実空間に落とせ、任意の $k$ に内挿できる。
戻すための $A(k)=S^{\rm PMT}(k)\,z^{\rm MLO}(k)$ は**保存済みのもの**を使う。
これで**内挿できないものを内挿しない**という条件を満たす。

### 必要な k の網羅 — ここが現状との差

$z$ を作る場所は現在 `sugw` 段 a' のみで、**GW の q リスト**（既約 BZ ＋ offset-Γ）しか回らない。
読む側は 3 種類ある:

| 読む場面 | 必要な $k$ | 現状 |
|---|---|---|
| 次反復の GW ドライバ | GW の q リスト | **覆えている** |
| 同反復の最後の `lmf`（SCF）| SCF の k メッシュ | `nkabc = n1n2n3` なら一致。**一般には一致しない** |
| バンドプロット | 経路上の任意の k | **覆えていない** → その場で作り直し ＝ 基底混在 |

**必要な変更は、$z$ を作る k を a' の q リストから「後で読む全ての k」へ広げること。**
$\Sigma_n$ を作るのと同じ `lmf --jobgw=1` の中で

1. GW の q リスト（既存）
2. SCF の k メッシュ
3. **バンド経路の k**（`syml` を先に決めておけば列挙できる）

の 3 つすべてで $z$ を計算し、`ZmloSig` に入れてしまう。
バンド経路を事前に固定する必要はあるが、診断用の経路なら問題にならない。

### これで何が保証されるか

- $\Sigma_n$ は**常に、それを作った $H$ の $\tilde\chi$** で読み戻される（任意の $k$ で）。
- $\tilde\chi$ は**反復ごとに更新される**（窓が古いバンド位置に張り付かない）。
- 内挿するのは $\Sigma^{\rm MLO}(R)$ だけ。$\langle\chi^{\rm PMT}|\tilde\chi\rangle$ は内挿しない。

### 残る近似

反復 $n$ の最後の `lmf` は SCF で密度を動かすので、その中で $H$ は変わる。
$z$ は反復開始時の $H$ のものを使い続けるため、SCF ループ内では厳密には整合しない。
ただし
- 同じ $z$ を使い続けるので **$\Sigma$ の読み方は一貫している**（揺れない）
- 自己無撞着に達すれば密度が止まるので、この近似は消える

つまり **途中の反復でのみ残る近似**であり、収束を妨げる種類のものではない、というのが見立て。

### 実装量

- `sugw` 段 a' を「q リスト」から「必要な k 集合」へ拡張（k の列挙と、そこでの `Hreduction` 呼び出し）
- SCF の k メッシュとバンド経路を `lmf --jobgw=1` の時点で知る必要がある
  （前者は `nkabc` から、後者は `syml.<sname>` から）
- 少し面倒だが、機構自体（`ZmloSig`/`ZmloNew` と昇格）は既にある

---

# 5. 過去の失敗との違い — 4 段階ある

| | 挙動 | 結果 |
|---|---|---|
| (a) 当初 | `getsenex` の**呼び出しごと**に作り直し。SCF ステップごとに基底が揺れる。読み書きの一致という概念すら無い | iter 2 で破綻 |
| (b) 凍結 `01b9ea4df` | 永久固定。しかも `getsenex` だけが追記するので iter 1（SCF メッシュ）と iter 2（GW の q）の**継ぎ接ぎ** | iter 1・2 は正常化、iter 3 で破綻 |
| (c) `5885ea2bb` | `lmf` 1 回の中で固定、次の `lmf` で作り直し。churn と継ぎ接ぎは消えた | **読み基底が 1 反復ずれたまま** |
| (d) **`b3dfecc89`（現行）** | 上記 2 スロット。$\tilde\chi$ は毎反復更新、$\Sigma$ は書かれた基底でのみ読む | **未検証** |

# 6. 残っている論点

**MLO 索引の窓（`eferm`, `ecbot`, `nskip`）は 0d で固定されたまま**です。
$z$ を今の $H$ で作り直しても、窓が初回のバンド位置を指していれば部分空間は古いままです。
ただし `ix` と `ndimMTO` は `SigRsMLO` との整合のため**変えられません**。
窓だけ現在の計算に追随させるかどうかは未決です。

---

# 7. 図

## 7.1 いまの比較（従来 MTO vs MLO 2 スロット版）

[![mlo_rows](Samples/kBT/LiTi2O4/mlo_rows.png)](Samples/kBT/LiTi2O4/mlo_rows.png)

（生成: [mlo_rows.py](Samples/kBT/LiTi2O4/mlo_rows.py)。左＝従来 pwmode=1、右＝MLO pwmode=11、
行は MLO の反復数で下に伸びる）

*表 7.1* t2g b33–44、Γ→X 211 点

| | 帯幅/本 | 弦ずれ | 弦/帯幅 | 前段からの変化 |
|---|---|---|---|---|
| LDA（内挿なし＝下限）| 391.6 meV | 30.8 | 7.9 % | — |
| 従来 MTO iter 1 | 290.7 | 28.7 | 9.9 % | 194.4 |
| 従来 MTO iter 2 | 223.0 | 25.4 | 11.4 % | 86.1 |
| **MLO 2 スロット iter 1** | 234.7 | 43.3 | **18.4 %** | 414.7 |
| **MLO 2 スロット iter 2** | 370.8 | 86.9 | **23.4 %** | 157.8 |

$E_{\rm HF}$ の iter 1→2 の動き: 従来 **2.85 eV**、MLO **8.73 eV**。

**注意**: この 弦/帯幅 には §4.2 ⑦ の**基底混在**が混ざっている。
メッシュ点では `ZmloSig` の $z$、その間はその場の $H$ から作った $z$ を使うため、
k 方向に基底の出どころが切り替わる。切り分けは `band3.sh`（混在／全 k 作り直し／
GW ドライバ状態の 3 通り）で行う。

## 7.2 バグの効果を分離した図（iter 1）

[![frozen_vs_not](Samples/kBT/LiTi2O4/frozen_vs_not.png)](Samples/kBT/LiTi2O4/frozen_vs_not.png)

（$\tilde\chi$ 凍結の有無。凍結なしだと帯幅が 292 → 202 meV に潰れる）

## 7.3 $\Sigma(R)$ の減衰（オンサイト規格化）

[![sigr_decay](Samples/kBT/LiTi2O4/sigr_decay.png)](Samples/kBT/LiTi2O4/sigr_decay.png)

（絶対値は正規化が違うので比較不可。規格化すると**従来 $\Sigma^{\rm MTO}(R)$ の方が減衰しない**）

---

# 8. 新方式 — 実装仕様

§4.6 の結論を実装できる形にしたもの。**これを実装する。**

## 8.1 不変条件

| | 内容 |
|---|---|
| **I1** | $\Sigma^{\rm MLO}_n$ は**必ず $z_n$ で読み戻す**。読む $k$ が何であっても例外を作らない |
| **I2** | $z_n$ は**$\Sigma_n$ を作ったハミルトニアン**から作る（反復 $n$ 開始時点の $H$）|
| **I3** | **内挿するのは $\Sigma^{\rm MLO}(R)$ だけ**。$\langle\chi^{\rm PMT}\vert \tilde\chi\rangle$ は内挿しない（APW 行が q 周期的でないため不可能）|
| **I4** | MLO **索引**（`ix`, `ndimMTO`）は連鎖を通じて固定。$\Sigma^{\rm MLO}$ がそのチャネル上の行列だから |

I1 と I2 を両立させるために、$z$ を**後で読む全ての $k$** で先に作っておく。これが新方式の核心。

## 8.2 k 集合 $K$

$$K = K_{\rm GW}\ \cup\ K_{\rm SCF}\ \cup\ K_{\rm band}$$

| | 中身 | 出どころ |
|---|---|---|
| $K_{\rm GW}$ | GW が要求する q（既約 BZ ＋ offset-Γ）| `QGpsi` / `m_qplist` |
| $K_{\rm SCF}$ | SCF の k メッシュ | `nkabc`（`n1n2n3` と一致するとは限らない）|
| $K_{\rm band}$ | バンド経路上の k | `syml.<sname>` から生成 |

**$K$ は `lmf --jobgw=1` の中で列挙できなければならない。**
$K_{\rm band}$ のためにバンド経路を事前に固定する必要があるが、診断用の固定経路なので問題ない。
`syml` が無ければ $K_{\rm band}=\emptyset$ とし、バンドは描かない（描くなら `syml` を置く）。

### 8.2.1 【後回し】MLO 空間で描けば $K_{\rm band}$ は不要になる

バンドを **PMT 空間**で描く限り、経路上の各 $k$ で $A(k)=S^{\rm PMT}(k)z^{\rm MLO}(k)$ が要るので
$K_{\rm band}$ を先に回しておく必要がある。

しかし **MLO 空間で描けばその必要が無い**。$H^{\rm MLO}(R)$ も $\Sigma^{\rm MLO}(R)$ も
**どちらも MLO 基底上の行列で q 周期的**なので、両方を任意の $k$ にフーリエ内挿して

$$\bigl[H^{\rm MLO}(k)+\Sigma^{\rm MLO}(k)\bigr]\,c = \varepsilon\,O^{\rm MLO}(k)\,c$$

を $n_{\rm MLO}\times n_{\rm MLO}$ で解けばよい。**PMT への引き上げが一切要らない。**

- $K = K_{\rm GW}\cup K_{\rm SCF}$ だけになる（$K_{\rm band}$ 不要）
- §8.5 の 3 点目（`syml` を `--jobgw=1` の前に用意）も不要になる
- 機構は既にある（`job_mlo`、`HamRsMLO` の `hammr`/`ovlmr`）

描かれるのは「縮約模型のバンド」であって PMT ハミルトニアンのバンドそのものではないが、
**MLO 模型が $\Sigma$ を再現できているかを見るのが目的**なので、診断としてはむしろ素直。

**2026-09-25 時点では後回し**（user 判断）。まず $K=K_{\rm GW}\cup K_{\rm SCF}$ で実装し、
バンドは従来どおり PMT 空間で描く（そのぶん経路上は作り直しになる）。

## 8.3 ファイル

| ファイル | 中身 | 寿命 |
|---|---|---|
| `SigRsMLO` | $\Sigma^{\rm MLO}(R)$ ＋ 対リスト ＋ `ix`, `ib_tableM` | 反復をまたぐ |
| `ZmloSig.<procid>` | **$k\in K$ 全部**での $z_n(k)$。現行 `SigRsMLO` と対 | 反復をまたぐ |
| `ZmloNew.<procid>` | 作成中の $z_{n+1}(k)$ | 反復内 |
| `HamRsMLO` 末尾 | `ix`, `ndimMTO`, `mlomethod`, `nskip`, `eferm`, `ecbot` | 連鎖を通じて固定 |

`ZmloSig` の被覆が $K_{\rm GW}$ だけから $K$ 全体に広がるのが、現行との唯一の構造差。

## 8.4 1 反復の手順

| # | 実行 | $\tilde\chi$ について |
|---|---|---|
| ① | `lmf --jobgw=1` | |
| | ├ `getsenex` | `ZmloSig` の **$z_{n-1}(k)$** で $\Sigma_{n-1}$ を読む。$k\in K$ なので**必ず見つかる** |
| | └ 段 a' | **$K$ の全ての $k$** でその場の $H$ から $z_n(k)$ を作り `ZmloNew` へ。$c^{\rm MLO}$ は $K_{\rm GW}$ でのみ計算して `__cmlo.data` へ |
| ② | GW 諸段 | $\Sigma^\psi$（MLO と無関係）|
| ③ | `hqpe_sc` | $\Sigma^{\rm MLO}_n=c_n^\dagger\Sigma^\psi c_n$ → `__SigmMLO.q`。**$z_n$ の基底** |
| ④ | `mlo --mlofreeze --mlo` | 対称化 → 全 BZ → FFT → `SigRsMLO`。**$z_n$ の基底** |
| ⑤ | **昇格** | `ZmloNew.*` → `ZmloSig.*` |
| ⑥ | `lmf`（Σ 入り SCF）| `ZmloSig` の $z_n$ で読む。$K_{\rm SCF}\subset K$ なので作り直し無し |
| ⑦ | `job_band ... --mlo` | 同じく $z_n$。$K_{\rm band}\subset K$ なので**作り直し無し**（基底混在が消える）|

### 8.4.1 ①の `getsenex` の式

既知なのは **MLO で挟んだ $\Sigma$**、すなわち
$\Sigma^{\rm MLO}_{\alpha\beta}=\langle\tilde\chi_\alpha|\hat\Sigma|\tilde\chi_\beta\rangle$ を実空間にしたもの。
これを PMT 基底の行列に戻す。**MLO は非直交なので $O^{-1}$ が両側に要る。**

**1. ブロッホ和**（唯一の内挿。$\Sigma^{\rm MLO}$ は MLO 基底上の行列なので実空間に落とせる）

$$\Sigma^{\rm MLO}_{\alpha\beta}(k)=\sum_{R}\Sigma^{\rm MLO}_{\alpha\beta}(R)\,e^{ikR}$$

**2. 引き上げ行列と重なり**（保存済みの $z^{\rm MLO}(k)$ ＋ その場の $S^{\rm PMT}(k)$）

$$A_{m\alpha}(k)=\langle\chi^{\rm PMT}_m|\tilde\chi_\alpha\rangle=\bigl[S^{\rm PMT}(k)\,z^{\rm MLO}(k)\bigr]_{m\alpha},
\qquad
O^{\rm MLO}_{\alpha\beta}(k)=\langle\tilde\chi_\alpha|\tilde\chi_\beta\rangle=\bigl[z^\dagger A\bigr]_{\alpha\beta}$$

**3. 非直交基底での射影演算子**

$$\hat P=\sum_{\alpha\beta}|\tilde\chi_\alpha\rangle\,(O^{-1})_{\alpha\beta}\,\langle\tilde\chi_\beta|
\ \Longrightarrow\
\hat P\hat\Sigma\hat P=\sum_{\alpha\beta\gamma\delta}|\tilde\chi_\alpha\rangle (O^{-1})_{\alpha\gamma}\,
\Sigma^{\rm MLO}_{\gamma\delta}\,(O^{-1})_{\delta\beta}\langle\tilde\chi_\beta|$$

**4. PMT 基底で挟む**

$$\mathrm{senex}_{mn}(k)=\langle\chi^{\rm PMT}_m|\hat P\hat\Sigma\hat P|\chi^{\rm PMT}_n\rangle
=\Bigl[A\,(O^{\rm MLO})^{-1}\,\Sigma^{\rm MLO}\,(O^{\rm MLO})^{-1}A^\dagger\Bigr]_{mn}$$

これが $H(k)$ に足される。設計書の式 (12)(14)(15) に対応。

**新方式で変わるのはここではない。** 式は同じで、**$z^{\rm MLO}(k)$ をどこから持ってくるか**だけが変わる
（その場で作り直す → $\Sigma$ を作った $H$ で作って保存したものを引く）。

## 8.5 現行コードとの差分

1. **`sugw` 段 a'** — ループを $K_{\rm GW}$ から $K$ 全体へ。
   $K_{\rm SCF}$ は `nkabc` から、$K_{\rm band}$ は `syml` から列挙する。
   $c^{\rm MLO}$（`__cmlo.data`）は従来どおり $K_{\rm GW}$ のみでよい。
2. **`m_sigmlo`** — $k$ が `ZmloSig` に無かったとき、いまは**黙って作り直す**。
   新方式では**それが起きたら異常**なので、**警告して止める**（`ECALJ_MLO_ALLOW_REBUILD=1` で従来動作）。
   混在を静かに見逃さないための歯止め。
3. **`gwsc`** — `syml.<sname>` を `lmf --jobgw=1` の前に用意する。

## 8.6 検証（これが揃って初めて「新方式で測った」と言える）

| 見るもの | 期待 |
|---|---|
| `llmfgw01` / `llmf` / `llmf_band` の「作り直し」回数 | **0**（1 回でも出たら $K$ の列挙漏れ）|
| `gwsc: promoted N ZmloNew.* -> ZmloSig.*` | 毎反復 |
| `loaded ZmloSig ... records=` | 各実行で $\vert K\vert $ 相当 |
| `ECALJ_ZMLO_DUMP` | a' と `getsenex` の $z$ が一致 |
| `ECALJ_SIGMLO_RT` | $\Sigma^{\rm MLO}(q)\to R\to q$ の往復誤差 |

## 8.7 残る近似

反復 $n$ の最後の `lmf` は SCF で密度を動かすが、$z$ は反復開始時のものを使い続ける。
したがって SCF ループ内では $z$ と $H$ が厳密には対応しない。ただし

- 同じ $z$ を使い続けるので **$\Sigma$ の読み方は一貫している**（k でも SCF ステップでも揺れない）
- 自己無撞着に達すれば密度が止まるので、この近似は消える

**途中の反復にのみ残る近似**であり、収束を妨げる種類のものではない、というのが見立て。

## 8.8 コスト

$K_{\rm band}$ ぶんの `Hreduction` が増える。1 k あたり PMT の対角化 1 回（$O(n_{\rm dimh}^3)$）。
211 点の経路なら GW ドライバの実行時間に数分乗る程度で、GW 諸段に比べれば小さい。

---

# 9. Q&A

## Q1. §8.4.1 の $O^{-1}$ は「部分空間の逆」か「PMT 空間の逆」か。どちらがよいか

**候補は 2 つ。**

| | 形 | 次元 |
|---|---|---|
| **(a) 部分空間の逆**（現実装）| $O^{\rm MLO}=\langle\tilde\chi\vert \tilde\chi\rangle=z^\dagger S z$、$\hat P=\vert \tilde\chi\rangle (O^{\rm MLO})^{-1}\langle\tilde\chi\vert $ | $n_{\rm MLO}\times n_{\rm MLO}$ |
| **(b) PMT 空間の逆** | $S^{-1}$ を使う形 | $n_{\rm dimh}\times n_{\rm dimh}$ |

**MLO については (a) しかない。** 理由は 2 つ。

1. $\tilde\chi_\alpha$ は **PMT 基底関数そのものではなく線型結合**である。
   従来の MTO 経路では部分空間が PMT の**部分ブロック**なので「$S$ の小行列を逆にする」ことに
   意味があり、実際 [getsenex](SRC/subroutines/rdsigm2.f90#L67) は `ovlm(1:ndimsig,1:ndimsig)` を
   逆にしている。**MLO には逆にすべき小行列が存在しない。**
2. $S$ を計量とする非直交基底で $\hat P^2=\hat P$ を満たす射影子は
   $|\tilde\chi\rangle O^{-1}\langle\tilde\chi|$ に**一意に決まる**。
   $S^{-1}$ を使う形は、$\Sigma^{\rm MLO}$（$n_{\rm MLO}^2$ 個の拘束）から
   PMT 係数 $\sigma$（$n_{\rm dimh}^2$ 個）を決めることになり**劣決定**で、
   最小性条件を課すと結局 (a) に戻る。

**未確認の点（正直に）**: 従来経路で `hqpe_sc` が書く `sene` の規約。
[main_hqpe.sc.f90](SRC/subroutines/main_hqpe.sc.f90) は**擬似逆**
$\psi^+=(\psi^\dagger\psi)^{-1}\psi^\dagger$ で挟み、`getsenex` 側は **$S_{\rm sub}^{-1}$** で挟み返す。
この 2 つが整合して $\langle\chi|\hat\Sigma|\chi\rangle$ を再現するのか、代数的に追い切れていない
（オンサイトが 73.8 eV という大きさは、行列要素そのものではなく双対側の量であることを示唆する）。

**数値で決着できる**: `ECALJ_SIGMLO_CHECK=1` は各 q で MLO 版と MTO 版の `senex` を直接比較する。
規約が食い違っていれば**系統的な因子**として現れる。
2026-09-25 に走らせたときは「全 q で一様に 13〜23 %」だったが、
あれは基底不整合のあるバイナリでの測定なので、**新方式で測り直す**（§8.6 の検証項目に追加）。

## Q2. なぜ従来法にはこの面倒が無いのか

引き上げ行列 $A$ が**その場の重なり行列の切り出し**だから。
従来の部分空間は PMT 基底の部分ブロックそのものなので $A=S^{\rm (full,sub)}$ であり、
保存すべきものが何も無い。MLO は $A=S^{\rm PMT}(k)\,z^{\rm MLO}(k)$ で、
$z$ が **$k$ 依存かつ $H$ 依存**である。ここが構造的な非対称。

## Q3. $z^{\rm MLO}$ を実空間に落として内挿すればよいのでは

**できない。** $z^{\rm MLO}$ の入った行列は **PMT 基底の行**を持ち、
PMT のうち q 周期的なのは MTO 部分だけで **APW 部分は q 周期的でない**
（`pwmode=11` は $|q+G|$ カットなので G の集合自体が q で変わる）。
行が $k$ ごとに別の関数を指すので実空間に落とせない。
同じ理由で $H^{\rm PMT}(R)$ を保存して $H^{\rm PMT}(k)$ を再構成する道も無く、
**$H^{\rm PMT}(k)$ はその密度で実際にハミルトニアンを組む以外に入手できない**。
→ だから §8 の方式（読む $k$ を全部先に回しておく）になる。

## Q4. 1 反復で MLO を何回作っているのか

**3 箇所**（2026-09-25 現在）。

| 場所 | メッシュ上の $q$ | メッシュ外の $k$ |
|---|---|---|
| `getsenex`（読む側）| `ZmloSig` から再利用 | **その場の $H(k)$ から作り直し** ← 新方式で禁止する |
| `sugw` 段 a'（書く側）| 毎回作り直し（意図的）| — |
| `mlo --mlofreeze --mlo` | 毎回作り直し。しかも **`__HamiltonianPMT` = 段 0c の LDA の $H$** から | — |

3 番目は `--mlofreeze` で `HamRsMLO` を書き直さないため、
**連鎖の最初の LDA ハミルトニアン**から毎回作り直している。
$\Sigma$ の対称化・回転・FFT は索引レベルの操作なので $\Sigma$ 自体は汚していないと思われるが、
**未確認**。新方式では確認する。

## Q5. 窓 $\bar\theta$ を広げれば $\Sigma(R)$ は短距離になるか

**ならない。** 2026-09-25 に「窓が狭いのでは」と書いたが撤回した。

- `mlomethod=4` では第 1 項が恒等的に 0（[`efrz = -1d99`](SRC/subroutines/m_hreduction.f90#L246)、
  「active only for mlomethod=3」）。効くのは第 2 項だけ。
- 第 2 項の位置は $\varepsilon^{\rm cut}_j=\max(\varepsilon_{\rm cbot}+\Delta,\ \varepsilon^{\rm MTO}_j)$。
  MTO のみの固有値は PMT のそれより上にあるので通常は $\varepsilon^{\rm MTO}_j$ が選ばれ、
  引数は $(\varepsilon^{\rm PMT}_i-\varepsilon^{\rm MTO}_j)/w$ という**バンド差**になって $k$ について滑らか。
- 幅も狭くない。LiTi₂O₄ は `mlo_w = 2.0` eV、`mlo_delta = 2.0` eV。

$\hat P(k)$ の $k$ 依存性は窓ではなく、
**固有ベクトル $|\Psi^{\rm PMT}_i(k)\rangle\langle\Psi^{\rm PMT}_i(k)|$ 自体が $k$ で回ること**から来る。

## Q6. $\tilde\chi$ にゲージ（位相）の任意性は入らないか

**入らない。** [m_hreduction.f90:322](SRC/subroutines/m_hreduction.f90#L322) で
$|F^{\rm MLO}_k\rangle=\hat P\,|F^{\rm MTO}_k\rangle$ の形になっており、
$\Psi^{\rm MTO}_j$ も $\Psi^{\rm PMT}_i$ も $|\Psi\rangle\langle\Psi|$ の形でしか現れないので
対角化の任意位相は相殺する。種は**裸の MTO 基底関数** $\chi^{\rm MTO}_k$（$k$ 非依存の実空間関数の
ブロッホ和）。
