# MLO-QSGW の理論確認（2026-09-25）

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
| ① | **`lmf --jobgw=1`**（`sugw`） | |
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
渡します（`SRC/exec/job_band`:91,93）。`--mlo` が無いと `c0_mlo` が偽になり、`sigmlo_init` は
冒頭の `if(.not.c0_mlo) return` で即座に戻ります。**MLO 経路は完全に無効**になり、`getsenex` は
従来の `sigm` を使う。`hqpe_sc` は `sigm` と `__SigmMLO.q` を並行して書くので `sigm` は常に
存在し、**黙って従来内挿のバンドが描かれます**。今日の図はすべてこれでした。

**(2) 連鎖ディレクトリの中で走らせない。** `job_band` の第 1 段は `lmf --quit=band` で、
これは **`efermi.lmf` を書き換えます**（`m_bndfp.f90`:274）。MLO の窓は `efermi.lmf` から
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
退避して a' で使います（`sugw.f90`）。2 スロット化しているので、`getsenex` 側は `ZmloSig` を
読むだけであり、a' がどの $H$ を使っても読み戻しは壊れません。

**これは標準の MLO 手法と同じです。** `m_bandcal.f90` では

| 行 | 処理 |
|---|---|
| 207 | `hambl` → $H^{\rm LDA}$ |
| 210–211 | `sigmamode` なら `getsenex` → **`hamm = hamm + senex`** |
| 221〜 | `if(writeham)` → **`__HamiltonianPMT` に書く** |

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
| `ZmloSig.<procid>` | 現行 `SigRsMLO` が書かれた $\tilde\chi$。`getsenex` が読む。追記しない |
| `ZmloNew.<procid>` | 段 a' がこの反復の $H$ から作った $\tilde\chi$。⑤で `ZmloSig.*` に rename |

ランクごとのファイルなのは、段 a' が MPI で $q$ を分担するのと、
`pwmode=11` で `ndimh` が $q$ ごとに変わり固定レコード長が使えないためです。

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
