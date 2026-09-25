# 自己エネルギーの MLO 表現による内挿 — 設計書

2026-09-24 起草、2026-09-25 に **§9〜§13 を追加**（`temp.md` を解体して統合）。対象: ecalj 開発者。
背景データは [kBT_research.md](kBT_research.md) の各日エントリ（図表番号はそちらを参照）。

> **現行の方式は §9〜§13**（2 スロット方式）。§3「設計」と §4「実装手順」は 2026-09-24 版で、
> 一部が置き換わっている。どこが置き換わったかは §9 冒頭の表にまとめてある。
> 実装の TODO は §13。

---

## 1. 何を解決するのか

QSGW のバンドに、**Σ の q メッシュ点の「間」にだけ現れる凸凹**が乗る。LiTi₂O₄ の E_F 直上
1 eV 付近で 25〜29 meV、O 2p 帯でも 13 meV 級。以下は実測で確定した事実である。

| 事実 | 根拠 |
|---|---|
| Σ メッシュ点の**上**では全バンドが滑らか。山は必ず隣接メッシュ点の**中点**に立つ | 図 22:00-1、14:00-1 |
| メッシュを密にしても振幅は落ちない（占有 t2g はむしろ増える） | 表 16:50-1 |
| `nkabc` を GW メッシュに一致させても**バンドは変わらない**（rms 1.5 meV） | 表 20:49-5、図 20:49-3 |
| 補間は既に Wigner–Seitz 最短ベクトル方式（局在を尊重）である | `m_shortn3.f90:gennlat`、`rdsigm2.f90:bloch2` |
| Σ(R) は Ti 3d ブロックだけが BvK セル端で頭打ち（0.3〜0.6 meV で減らない） | 09-24 の Σ(R) 測定 |
| 生の MTO 表現で自由度を削ると Σ の値が数百 meV 動く／発散する | `ECALJ_SIG_1RAD`, `ECALJ_SIG_EHONLY` の試験 |
| MLO 表現なら**値を保ったまま**自由度を減らせる（メッシュ点で 0.1〜0.5 meV） | 表 00:05-1、表 22:45-1 |

原因は**自由度過多**である。6³ なら既約 q は 16 点しかないのに、実空間では 216 セル × 230² の
自由度を持たせている。しかも MTO 基底は非直交で EH と EH2 が準線形従属、裾も長い。

MLO 表現（154 次元、局在、固定種でゲージ固定）に移すと、Γ–X の補間の荒れは
**9.8 → 6.3 meV**（全 EH lm 154 軌道）、**9.8 → 3.5 meV**（t2g 12 軌道）に落ちる。

---

## 2. 現状のデータフロー — 出発点から順に

### 2.1 GW が出すもの

`hsfp0` / `hgw` が各既約 $q$ で計算するのは、**固有関数で挟んだ行列要素**である:

$$\Sigma^{\psi}_{ij}(q) \;=\; \langle \psi^{\rm PMT}_{iq} \,|\, \hat\Sigma \,|\, \psi^{\rm PMT}_{jq}\rangle\tag{1}
$$

$i,j$ はバンド指標（`emax_sigm` 以下）。QSGW ではこれをエルミート化した静的 $\Sigma$ を使う。
固有関数は正規直交（$\langle\psi_i\mid\psi_j\rangle = \delta_{ij}$）なので、演算子は逆行列なしに

$$\hat\Sigma(q) \;=\; \sum_{ij} |\psi^{\rm PMT}_{iq}\rangle\;\Sigma^{\psi}_{ij}(q)\;\langle \psi^{\rm PMT}_{jq}|\tag{2}
$$

と組める。**ここが唯一の第一原理的な出発点**で、以下はすべてその表現の変換にすぎない。

### 2.2 `hqpe_sc`: MTO 基底の行列へ

固有関数を PMT 基底で展開する係数を $z_{mi}$ とすると
$|\psi_{iq}\rangle = \sum_m |\chi^{\rm PMT}_{mq}\rangle z_{mi}$、したがって
$\langle\chi^{\rm MTO}_\mu\mid\psi_i\rangle = \sum_m S_{\mu m}\,z_{mi}$
（$S_{\mu m} = \langle\chi^{\rm MTO}_\mu\mid\chi^{\rm PMT}_m\rangle$）。これで

$$\Sigma^{\rm MTO}_{\mu\nu}(q) \;=\; \langle\chi^{\rm MTO}_\mu|\hat\Sigma(q)|\chi^{\rm MTO}_\nu\rangle
\;=\; \sum_{ij}\Big(\sum_m S_{\mu m} z_{mi}\Big)\,\Sigma^{\psi}_{ij}(q)\,
\Big(\sum_n S_{\nu n} z_{nj}\Big)^{*}\tag{3}
$$

これが `sigm` の中身（$L\times L$、$L$ = `ldim` = 230）。**途中産物**である。

### 2.3 全 BZ へ、そして実空間へ

`rdsigm2` が対称操作で既約 $q$ から全 BZ の $q$ へ回転し（`hamfb3k`）、`fftz3` が実空間へ移す:

$$\Sigma^{\rm MTO}_{\mu\nu}(R) \;=\; \frac{1}{N_q}\sum_{q} \Sigma^{\rm MTO}_{\mu\nu}(q)\,e^{-iqR}\tag{4}
$$

### 2.4 任意 $k$ へ — ここが内挿

`bloch2` が Wigner–Seitz 最短ベクトル（`m_shortn3.f90:gennlat`、縮退した像は `nqwgt` で重み平均）で Bloch 和する:

$$\Sigma^{\rm MTO}_{\mu\nu}(k) \;=\; \sum_{R} \Sigma^{\rm MTO}_{\mu\nu}(R)\,\overline{e^{ikR}}\tag{5}
$$

**$\mu\nu$ 成分ごとに独立にフーリエ内挿される。**

### 2.5 `getsenex`: PMT 基底へ戻して $H$ に加える

MTO 基底は非直交なので、演算子は $\big(O^{\rm MTO}\big)^{-1}$ を両側に挟んだ形になる
（$O^{\rm MTO}_{\mu\nu} = \langle\chi^{\rm MTO}_\mu\mid\chi^{\rm MTO}_\nu\rangle$）:

$$\hat\Sigma \;=\; \sum_{\mu\rho\sigma\nu} |\chi^{\rm MTO}_\mu\rangle\,
\big(O^{\rm MTO}\big)^{-1}_{\mu\rho}\;\Sigma^{\rm MTO}_{\rho\sigma}(k)\;
\big(O^{\rm MTO}\big)^{-1}_{\sigma\nu}\,\langle\chi^{\rm MTO}_\nu|\tag{6}
$$

その PMT 基底での行列要素が `senex`:

$$\big[\hat\Sigma\big]_{mn}(k) \;=\; \sum_{\mu\rho\sigma\nu}
S^{*}_{\mu m}\,\big(O^{\rm MTO}\big)^{-1}_{\mu\rho}\;\Sigma^{\rm MTO}_{\rho\sigma}(k)\;
\big(O^{\rm MTO}\big)^{-1}_{\sigma\nu}\,S_{\nu n}\tag{7}
$$

コードの `ovliovl` が $\big[(O^{\rm MTO})^{-1}S\big]_{\mu n}$ に当たる。
`getsenex(qp,isp,ndimh,ovlm)` の呼び出しは **2 箇所だけ**（`m_bandcal.f90`:198,210 と `sugw.f90`:429,440）で、
どちらも引数 `ovlm` に $S^{\rm PMT}(q)$ を渡している。
$m,n$ は **MTO と APW の両方**を走るので **APW ブロックも埋まる**。
`mtosigmaonly` は「$\Sigma$ の行列要素を MTO 部分空間で保持する」の意であって、
「APW に効かない」ではない。

### 2.6 この経路が何をしているか

2.2 と 2.5 を合わせると、メッシュ点上でも

$$\hat\Sigma \;\to\; \hat P\,\hat\Sigma\,\hat P,
\qquad \hat P = \sum_{\mu\nu} |\chi^{\rm MTO}_\mu\rangle\big(O^{\rm MTO}\big)^{-1}_{\mu\nu}\langle\chi^{\rm MTO}_\nu|\tag{8}
$$

すなわち **MTO 部分空間への射影**が既に入っている（$\hat P$ は非直交基底に対する射影子）。
固有関数が MTO 部分空間に完全に含まれない限り、これは近似である。

**そして問題はここ**。6³ なら既約 $q$ は 16 点しかないのに、2.3–2.4 は
$L^2 = 230^2$ 個の要素それぞれについて $n_k^3 = 216$ セル分の実空間自由度を持たせている。
しかも MTO 基底は非直交で EH/EH2 が準線形従属、包絡関数の裾も長い。
結果として、**メッシュ点の間で線形独立性の低い自由度がふらつく**（§1 の表）。

```
hsfp0/hgw  ->  Sigma^psi_ij(q)        固有関数基底、既約 q
hqpe_sc    ->  Sigma^MTO_munu(q)      sigm
rdsigm2    ->  全 BZ (hamfb3k)
fftz3      ->  Sigma^MTO_munu(R)      hrr          <- 自由度 L^2 x nk^3
bloch2(k)  ->  Sigma^MTO_munu(k)                   <- 内挿はここだけ
getsenex   ->  [Sigma]_mn(k)          senex  -> H
```

## 3. 設計

### 3.0 記号 — 基底・MLO・固有関数を区別する

| 記号 | 意味 | 次元 / 索引 |
|---|---|---|
| $\chi^{\rm PMT}_m$ | **PMT 基底関数**（MTO + APW） | $m = 1\ldots n_{\rm dimh}$。APW の G の選び方は §3.2.1 の脚注 |
| $\chi^{\rm MTO}_\mu$ | その MTO 部分 | $\mu = 1\ldots L$（`ldim` = 230）、**$k$ 非依存** |
| $\chi^{\rm MTO}_{ix(\alpha)}$ | **MLO の種**。`ix(`$\alpha$`)` 番目の MTO 基底関数（固定） | 新しい記号は要らない |
| $\tilde\chi_\alpha$ | **MLO**（コードの `F^MLO`）。種をエネルギー窓付きで band 多様体へ射影したもの（式 (9)）。**冪等な射影ではない** | $\alpha = 1\ldots M$（`ndimMTO` = 154） |
| $\psi^{\rm PMT}_i$ | **固有関数**（基底ではない） | |
| `cmlo`$_{i\alpha} = \langle \psi^{\rm PMT}_i\mid\tilde\chi_\alpha\rangle$ | MLO の**固有関数**基底での係数 | $(n_{\rm band}, M)$、**両添字とも $q$ 非依存の数** |
| `zMLO`$_{m\alpha}$ | MLO の **PMT 基底**での係数、$\tilde\chi_\alpha=\sum_m \chi^{\rm PMT}_m z^{\rm MLO}_{m\alpha}$ | $(n_{\rm dimh}, M)$ |

コードの `F^MTO_k` は $\chi^{\rm MTO}_{ix(\alpha)}$、`F^MLO_k` は $\tilde\chi_\alpha$ に対応する。
**MLO の索引には $\alpha,\beta$ を使い、波数 $k,q$ と区別する。**

$\Sigma$ は `sigm` に $\Sigma^{\rm MTO}_{\mu\nu}(q)=\langle\chi^{\rm MTO}_\mu|\hat\Sigma|\chi^{\rm MTO}_\nu\rangle$ として入っている（$L\times L$）。

### 3.1 変えられるのは部分空間だけ

§2.5–2.6 のとおり、展開の枠組み（$\big(O\big)^{-1}$ を両側に挟む非直交の射影）は既に一般的な形であり、
**APW ブロックも埋まっている**。したがって設計の自由度は

> **どの部分空間で $\Sigma$ を保持し、内挿するか**

の一点に尽きる。現状は MTO（$L = 230$、非直交、EH/EH2 が準線形従属、裾が長い）。
答は **MLO**（§3.2）。式 (1) の $\Sigma^{\psi}$ から直接その部分空間の行列を作り、**MTO を経由しない**。

なお **MT 球の内外への分割は、この問題に関係しない**。球内の部分波チャネルで $\Sigma$ を保持する道も考えたが、
チャネル数は `ndima` = 700 で MTO の 230 より**増える**うえ、band 多様体の窓で絞る仕組みが無いので、
自由度過多という §1 の病そのものに効かない。格子間の $\Sigma$ も落ちる。採らない。

### 3.2 MLO を部分空間に使う

#### 3.2.1 MLO は PMT 基底の線型結合である

本設計が寄りかかっている事実はこれ一つである。**MLO は PMT 基底関数の線型結合として書かれており、
その係数はコードが既に作っている**（`Hreduction` の `zMLO`、`cmlo`）:

$$\tilde\chi_{\alpha q} \;=\; \sum_m \chi^{\rm PMT}_{mq}\, z^{\rm MLO}_{m\alpha}(q)
\;=\; \sum_i \psi^{\rm PMT}_{iq}\, c^{\rm MLO}_{i\alpha}(q),
\qquad z^{\rm MLO} = z^{\psi}\, c^{\rm MLO}\tag{9}
$$

$z^{\psi}_{mi}$ は固有ベクトル（`evecpmt`）。$\psi$ は正規直交なので
$c^{\rm MLO}_{i\alpha}=\langle\psi^{\rm PMT}_{iq}\mid\tilde\chi_{\alpha q}\rangle$ そのものである。

**ブロッホ基底とは「$q$ について周期性を持つもの」である。**
すなわち $\chi_{\mu,q+G} = \chi_{\mu q}$（$G$ は逆格子ベクトル）。
これは「**$q$ に依らない実空間関数 $f$ ひとつからブロッホ和 $\sum_T e^{iqT} f(\mathbf r-\mathbf T)$ で構成できる**」
ことと同値である。片方向は $e^{i(q+G)T}=e^{iqT}$ から直ちに出る。逆も言えて、$q$ の周期関数なら
$f(\mathbf r-\mathbf T)=\frac{1}{N_q}\sum_q e^{-iqT}\chi_{\mu q}(\mathbf r)$ で種が取り出せ、和を取れば元に戻る。
**これが $q$ のフーリエ級数＝実空間表現が存在する根拠そのもの**である。

$$\chi^{\rm MTO}_{\mu q}(\mathbf r) = \sum_{T} e^{iqT}\,\chi^{\rm MTO}_{\mu}(\mathbf r - \mathbf T - \boldsymbol\tau_{ib})
\quad\Rightarrow\quad \chi^{\rm MTO}_{\mu,q+G} = \chi^{\rm MTO}_{\mu q}\tag{10}
$$

- **$\chi^{\rm MTO}$ はブロッホ基底**。だから $\Sigma^{\rm MTO}_{\mu\nu}(q)$ は $q$ の周期関数で、
  式 (4)(5) のフーリエ変換が意味を持つ。`hammr`/`ovlmr` も同じ理由で成り立つ。
- **$\chi^{\rm APW}_{Gq}=e^{i(q+G)\mathbf r}$ はブロッホ基底ではない**。$q\to q+G'$ で、ラベル $G$ の関数は
  ラベル $G'+G$ の関数に移る。**ラベルがずれる**ので、$G$ で番号付けた基底関数は $q$ の周期関数にならない。
  これは**本数の問題ではなくラベルの問題**である。

  > **APW の $G$ の選び方（`m_igv2x.f90`）**: `pwgmax = sqrt(pwemax)` として、
  > `pwmode` の 10 の位が **0 なら $|G| <$ `pwgmax`**（`qqq = 0`。**本数 `napw` は $q$ に依らない**）、
  > **1 なら $|q+G| <$ `pwgmax`**（$q$ 依存）。LiTi₂O₄ の計算は `pwmode = 1` なので**前者**。
  > どちらでも上のラベルずれは起きるので、$\mathrm{FFT}[z^{\rm MLO}]$ が取れない結論は変わらない。
  > なお厳密には、$|G|$ カットだと張る空間自体も $q$ の周期関数にならない（$G$ 球が原点固定だから）。
  > $|q+G|$ カットならそこは周期的になる。

**MLO はどうか。$\tilde\chi_\alpha$ は関数としては $q$ の周期関数である** — 種 $\chi^{\rm MTO}_{ix(\alpha)}$ が
ブロッホ基底、band 多様体も窓（エネルギーで決まる）も $q$ の周期的な量だから、$\tilde\chi_{\alpha,q+G}=\tilde\chi_{\alpha q}$。
**したがって $\Sigma^{\rm MLO}_{\alpha\beta}(q)$ は $q$ の周期関数であり、式 (11)(12) のフーリエ変換が正当化される。
これが本設計の licence である。**
同じ理由で $\tilde\chi_\alpha$ 自身も実空間の種を持ち、和を取れば任意 $k$ で再構成できる（式 (15) の $A(k)$ はそれ）。

周期的でないのは**係数配列 $z^{\rm MLO}_{m\alpha}(q)$ の APW 部分**の方で、$q\to q+G'$ でラベルがずれる。
だから $\mathrm{FFT}[z^{\rm MLO}]$ は取れない。しかし**内挿に要るのは行列であって係数ではない**（§3.2.3）。

| | |
|---|---|
| **構成法** | 固定した MTO 種 $\chi^{\rm MTO}_{ix(\alpha)}$ を、**エネルギー窓付きで** band 多様体へ射影する（`Amat` = 重なり × 窓の重み、`nskip` 以下は捨てる）。窓の具体形とパラメータ（`mlo_delta`, `mlo_wfrz`, `mlo_w`）は `m_hreduction.f90` を見ること。**冪等な射影ではない** |
| **窓** | $\Sigma^{\rm MLO}$ は**窓の中でのみ $\hat\Sigma$ を忠実に表す**。窓の外は鋭く切れるのではなく連続に落ちる。これは欠陥ではなく設計意図（§3.3） |
| **ゲージ** | 種が $q$ に依らず固定 ⇒ **射影ゲージ**。固有ベクトルの任意位相・縮退内の任意回転は相殺し、$z^{\rm MLO}(q)$ は $q$ の滑らかな関数（Wannier の gauge fixing が要らない） |
| **直交性** | 既定では**非直交**。$O^{\rm MLO}(q) = (z^{\rm MLO})^\dagger S^{\rm PMT}(q)\, z^{\rm MLO}$。`--mlo_ortho` のときだけ Löwdin 直交化 |
| **次元** | $\alpha = 1\ldots M$、$M$ = `ndimMTO`（LiTi2O4 全 EH lm で 154、t2g 模型で 12）。$i$ はバンド指標、$m$ は PMT 指標 |
| **出力** | 段 1 で `zMLO` の MTO 行を `__amlo.data` に書き出し済み。**この先で要るのは `cmlo` の方**（式 (11)）で、まだ外に出していない → 段 1' |

#### 3.2.2 式 (3)–(5) は MTO を MLO に置き換えるだけで済む

§2 の流れのうち **(3) ψ→基底 / (4) FFT / (5) Bloch 和** は、MTO を $\tilde\chi$ に
読み替えるだけで成立する。しかも式 (3) の MLO 版は **`zMLO` を経由しない**:

$$\Sigma^{\rm MLO}_{\alpha\beta}(q) \;=\; \langle \tilde\chi_\alpha|\hat\Sigma(q)|\tilde\chi_\beta\rangle
\;=\; \sum_{ij} \big(c^{\rm MLO}_{i\alpha}\big)^{*}\,\Sigma^{\psi}_{ij}(q)\,c^{\rm MLO}_{j\beta}
\;=\;\big((c^{\rm MLO})^\dagger\,\Sigma^{\psi}\,c^{\rm MLO}\big)_{\alpha\beta}\tag{11}
$$

$c^{\rm MLO}_{i\alpha}=\langle\psi^{\rm PMT}_{iq}\mid\tilde\chi_{\alpha q}\rangle$（式 (9)）の第一添字は
**バンド指標**なので、APW の事情は**ここには入らない**。
$\psi$ は正規直交だから逆行列も不要。$\Sigma^{\psi}$ から一発で取れる。

続く 2 段はそのまま:

$$\Sigma^{\rm MLO}_{\alpha\beta}(R) = \frac{1}{N_q}\sum_q \Sigma^{\rm MLO}_{\alpha\beta}(q)e^{-iqR},
\qquad
\Sigma^{\rm MLO}_{\alpha\beta}(k) = \sum_R \Sigma^{\rm MLO}_{\alpha\beta}(R)\,\overline{e^{ikR}}\tag{12}
$$

$\alpha,\beta$ は MLO の通し番号（$M$ = `ndimMTO`）で **$q$ 非依存**だから、FFT も Bloch 和も
MTO 版と同一の実装で回る。**内挿される自由度は $L^2 n_k^3 = 230^2\cdot n_k^3$ から
$M^2 n_k^3 = 154^2\cdot n_k^3$（t2g モデルなら $12^2\cdot n_k^3$）に減り、
しかも準線形従属な方向が取り除かれている** — §1 の表の問題に直接効く。

#### 3.2.3 実空間 MLO 表現は既にコードにある

$\tilde\chi$ は実空間表現を持ち、**それを使う機構は既に動いている**（MLO バンドプロットがそれ）。

| 実体 | 作る場所 | 中身 |
|---|---|---|
| `hammr`, **`ovlmr`** | `HamPMTtoHamRsMLO`（`m_HamPMT.f90`）が `HamRsMLO` に書く | $H^{\rm MLO}(R)$, $O^{\rm MLO}(R)$。添字は `npair(ib1,ib2)` / `nlat` / `nqwgt` の**対リスト**上。これは `bloch2` と同じ WS 最短ベクトル機構 |
| 任意 $q$ での Bloch 和 | `m_mlo_ham::calc_ham_eigen` | 上の 2 つを $q$ で和して $H^{\rm MLO}(q), O^{\rm MLO}(q)$ を作り一般化固有値問題を解く → `band_MLO_spin1.dat` |

$\Sigma^{\rm MLO}(R)$ は `hammr`/`ovlmr` と**まったく同種の MLO×MLO 行列**（添字 $\alpha\beta$ は $q$ 非依存）だから、
式 (12) は新規実装ですらなく、既存の経路に 1 本足すだけである。$O^{\rm MLO}(k)$ も既に任意 $k$ で得られている。

**APW の事情が効くのは 1 箇所だけ**: $z^{\rm MLO}_{m\alpha}(q)$ を**要素ごとにフーリエ内挿**しようとしたときである。
行索引 $m$ が PMT 基底なので $\mathrm{FFT}[z^{\rm MLO}]$ は定義できない。**その手を使わなければよい**だけで、
$z^{\rm MLO}$ の APW 部分が $q$ の周期関数でないからである（§3.2.1）。
**内挿するのは係数 $z^{\rm MLO}$ ではなく、$q$ の周期関数である行列 $\Sigma^{\rm MLO}_{\alpha\beta}$ の方**である。

#### 3.2.4 使い方は 2 通り

**(A) バンド / QP エネルギーが目的なら、追加実装はほぼ無い。**
$H^{\rm MLO}(R)$ に $\Sigma^{\rm MLO}(R)$ を足して `calc_ham_eigen` に渡すだけ。
$\Sigma$ を PMT 基底に戻す必要は無いので、式 (6)(7) に相当する段自体が存在しない。
09-24 の測定で、MLO 表現の Γ–X 粗さは MTO 内挿の 9.8 meV に対し 3.5 meV（t2g 12 軌道）/ 6.3 meV（154 軌道）、
メッシュ点の値は 0.1–0.5 meV で再現していた。

**(B) SCF（電子密度）まで回すなら、$\Sigma$ を PMT ハミルトニアンに戻す必要がある。**
lmf は PMT 基底で密度を作るので、式 (12) の $\Sigma^{\rm MLO}_{\alpha\beta}(k)$ を
ハミルトニアンの行列要素に書き換える式が要る。MLO も非直交なので
$\big(O^{\rm MLO}\big)^{-1}$ を両側に挟む（式 (6)(7) の MLO 版）:

$$\hat\Sigma \;=\; \sum_{\alpha\alpha'\beta'\beta} |\tilde\chi_{\alpha k}\rangle\,
\big(O^{\rm MLO}\big)^{-1}_{\alpha\alpha'}\;\Sigma^{\rm MLO}_{\alpha'\beta'}(k)\;
\big(O^{\rm MLO}\big)^{-1}_{\beta'\beta}\,\langle\tilde\chi_{\beta k}|\tag{13}
$$

$$\boxed{\;\big[\hat\Sigma\big]_{mn}(k) \;=\; \sum_{\alpha\alpha'\beta'\beta}
A_{m\alpha}(k)\,\big(O^{\rm MLO}\big)^{-1}_{\alpha\alpha'}\;\Sigma^{\rm MLO}_{\alpha'\beta'}(k)\;
\big(O^{\rm MLO}\big)^{-1}_{\beta'\beta}\,A^{*}_{n\beta}(k)\;}\tag{14}
$$

これが `senex(ndimh,ndimh)` である。要る行列は 2 つだけ:

$$A_{m\alpha}(k) = \big\langle\chi^{\rm PMT}_{mk}\mid\tilde\chi_{\alpha k}\big\rangle
= \big(S^{\rm PMT}(k)\,z^{\rm MLO}(k)\big)_{m\alpha},
\qquad
O^{\rm MLO}_{\alpha\beta}(k) = \big((z^{\rm MLO})^\dagger A\big)_{\alpha\beta}\tag{15}
$$

$m,n$ は MTO と APW の両方を走るので、現状（式 (7)）と同じく **APW ブロックも埋まる**。
$S^{\rm PMT}(k)$ は `getsenex` の引数 `ovlm` で既に渡っている。

**要件: 式 (15) は任意の $k$ で計算できなければならない。** `getsenex` の呼び出し元は 2 つあり、
**どちらでも**要る:

| 呼び出し元 | モード | $k$ リスト |
|---|---|---|
| `m_bandcal.f90` | 通常の lmf（密度・バンド） | `nkabc` の BZ メッシュ、バンドプロットの対称線 |
| `sugw.f90` | **lmf の GW ドライバモード**（`lmf --jobgw`） | GW の q リスト。既約 q に加え **offset-Γ 点 `q0p`** を含む |

後者を落とすことはできない。GW ドライバが書き出す固有関数・固有値が次反復の $\Sigma^{\psi}$ を決めるので、
ここで $\Sigma$ の入り方が通常モードと違えば、**反復そのものが不整合になる**。
逆に言えば、置き換えを `getsenex` の中で行う限り**両モードが自動的に同じ処方になる**のが、この設計の利点である。

どちらも `getsenex` は **$H$ を組み立てる途中**で呼ばれる。

これは満たせる。`getsenex(qp,isp,ndimh,ovlm)` が呼ばれる時点で

- `hamm` = $H^{\rm LDA}(k)$（`hambl` の直後。`senex` はこの呼び出しの**後**に足される）
- `ovlm` = $S^{\rm PMT}(k)$（引数そのもの）

がその場にある。これは `Hreduction(mlomethod,...,hamm,ovlm,...,qp,cmlo,nev,zMLO,...)` の入力そのものなので、
**その $k$ で `Hreduction` を呼べば $z^{\rm MLO}(k)$ が直接得られる**。$\mathrm{FFT}[z]$ は使わないので
APW の事情に関係なく、offset-Γ のような特殊点でも同じ手順で通る。

**循環しない理由**: MLO を $H^{\rm LDA}$ から作るから。$\Sigma$ を含む $H$ を使おうとすると
「$\Sigma$ を作るのに $\Sigma$ が要る」になるが、$H^{\rm LDA}$ なら `getsenex` 突入時点で確定している。
$\tilde\chi_\alpha$ は $\Sigma$ の**表現**であって $\Sigma$ の固有関数である必要は無いので、これでよい。

**そのために必要な変更（段 1' に反映）**: メッシュ点で $\Sigma^{\rm MLO}(q)$ を作るときの `cmlo` も
**同じく $H^{\rm LDA}$ から**作らねばならない。現状 `HreductionIqibz` は `sigm` があれば
$\Sigma$ 込みの $H$ で MLO を作るので、ここを LDA 側に固定する。
あわせて窓の大域量（`nskip`、$\varepsilon_{\rm cbot}$、$E_F$）は**メッシュ走査で決めた値を凍結**し、
`__cmlo.info` に書いて任意 $k$ でもそれを使う（$k$ ごとに決め直すと窓の意味が $k$ ごとに変わる）。

**コスト**: `Hreduction` は PMT ブロックと MTO ブロックを対角化するので、$k$ あたり対角化 2 回分が増える。

**では $A$ や $z^{\rm MLO}$ を内挿すれば安くならないか**（重い場合の逃げ道。**未検証**）。
$q$ 周期性（§3.2.1）で仕分けると、ブロックごとに事情が違う:

| 量 | $q$ 周期か | 内挿できるか |
|---|---|---|
| $\Sigma^{\rm MLO}_{\alpha\beta}$, $O^{\rm MLO}_{\alpha\beta}$ | **周期** | フーリエ内挿できる（本設計が使うのはこれ） |
| $A$ の **MTO 行** $\langle\chi^{\rm MTO}_{\mu k}\mid\tilde\chi_{\alpha k}\rangle$ | **周期**（両者ともブロッホ基底） | フーリエ内挿できる |
| $A$ の **APW 行** $\langle\chi^{\rm APW}_{Gk}\mid\tilde\chi_{\alpha k}\rangle$ | 周期でない（$G$ ラベルがずれる） | フーリエ内挿は**できない** |
| $z^{\rm MLO}$ の APW 行 | 同上 | 同上 |

ただし APW 行には別の構造がある。$\tilde\chi_{\alpha,k+G'}=\tilde\chi_{\alpha k}$ なので

$$A_{G\alpha}(k) \;=\; \tilde F_\alpha(\mathbf k+\mathbf G),
\qquad \tilde F_\alpha(\mathbf p) = \int d\mathbf r\; e^{-i\mathbf p\cdot\mathbf r}\,\tilde\chi_\alpha(\mathbf r)\tag{16}
$$

すなわち **$k$ と $G$ に別々に依存するのではなく、ベクトル $\mathbf k+\mathbf G$ ひとつの滑らかな関数**であり、
それは**実空間 MLO のフーリエ変換**に他ならない（$\tilde\chi_\alpha$ が局在しているから滑らかに減衰する）。
したがって「$k$ についてフーリエ内挿する」のではなく、
**$\tilde F_\alpha(\mathbf p)$ を逆空間の細かい格子に一度だけ作って参照する**のが筋になる。
球内の増補も $\mathbf k+\mathbf G$ で決まるので同じ扱いができるはず。**要検証**。

**整合条件**: 要するに $\tilde\chi_\alpha$ は、メッシュ点でも任意 $k$ でも**同一の処方**で作られた
同じ関数でなければならない（$H^{\rm LDA}$ + 凍結した窓）。
種 $\chi^{\rm MTO}_{ix(\alpha)}$ が固定なのでゲージは自動的に揃う（§3.2.1）。

### 3.3 何を失うか

捨てるのは **MLO 部分空間の外**、すなわち窓の外のバンド成分である。
LiTi₂O₄ で全 EH lm（154）を採れば窓内は 0.5 meV で再現されるが、
窓の外（E_F+4 eV 以上の格子間バンド）の $\Sigma$ は再現されない。DOS や光学応答を取るときは注意。

ただし、現状の式 (8) が既に MTO 部分空間への射影 $\hat P\hat\Sigma\hat P$ になっている以上、
「何かを捨てる」こと自体は新しくない。問題は**何を捨てるかを制御できるか**であり、
MTO では制御できない（§7 の付録）のに対し MLO では窓で制御できる、というのが本設計の要点である。

### 3.4 基底の $q$ 周期性と `pwmode` — 両立しないトレードオフ

§3.2.1 の licence は「$\Sigma^{\rm MLO}_{\alpha\beta}(q)$ が $q$ の周期関数であること」だった。
ところが **PMT 基底そのものの $q$ 周期性が `pwmode` の選択に依存する**。

| `pwmode` | G の選び方 | `napw` | 基底の $q$ 周期性 | $q$ についての連続性 |
|---|---|---|---|---|
| 1 | $\vert G\vert <$ `pwgmax` | $q$ 非依存 | **無い** | **滑らか**（基底が $q$ に依らない） |
| 11 | $\vert q+G\vert <$ `pwgmax` | $q$ 依存 | **有る** | 不連続（$\vert q+G\vert$ が cutoff を跨ぐたび基底が入れ替わる） |

**この 2 つは両立しない。** $q$ 非依存の有限集合 $S$ では
$E_n(q+G_i)$ が $\{G_i+S\}$ から作られるのに対し $E_n(q)$ は $\{S\}$ から作られ、
$G_i+S = S$ は有限集合では有り得ないので、**基底を $q$ 非依存に固定する限りバンドは BZ 周期的にならない**。
逆に周期性を取ると $|q+G|$ カットになり、$q$ を動かすと基底が飛ぶ。

実務上の使い分けもこの通りになっている:

- **QSGW（自己無撞着・$\Sigma$ の内挿）は $|q+G|$ カットであるべき**。
  基底が $q$ 周期でなければ $\tilde\chi_\alpha$ も $\Sigma^{\rm MLO}(q)$ も厳密には $q$ 周期でなく、
  式 (11)(12) のフーリエ変換に系統誤差が入る。§3.2.1 の licence が厳密に成り立つのは `pwmode = 11` のときだけ。
- **バンドプロットは $|G|$ カットでないと綺麗に描けない**。$|q+G|$ カットだと
  基底の入れ替わりが **BZ 内でジャンプ**として見える。基底を固定すれば $H(q), S(q)$ が $q$ の滑らかな関数になり、
  固有値も滑らかに動く。

**未解決（user 2026-09-24）**: $|G|$ カットでも、本当は $|G|$ と $|G+G_i|$（$G_i$ は隣接する逆格子ベクトル）の
**両方**で切らないと、バンド端で隣のゾーンと滑らかに繋がらないのではないか。
筋は通る — $S' = \{G : |G|<G_{\max} \text{ かつ } |G+G_i|<G_{\max}\}$ とすれば
$S'$ と $G_i+S'$ の食い違いが $G_{\max}$ 近傍（＝高エネルギー成分）だけに押しやられるので、
低いバンドのゾーン境界での一致は良くなるはず。ただし $G_i+S'=S'$ にはならないので**厳密な周期性は回復しない**。
**未検証**。

**本設計への含意**:
1. 段 2 で $\Sigma^{\rm MLO}(R)$ を作る計算は `pwmode = 11` で回すのが筋。
2. その $\Sigma^{\rm MLO}(R)$ を使ったバンドプロットは `pwmode = 1` でもよい
   （$\Sigma$ 側はもう実空間に落ちていて $q$ 周期だから、基底の選択と独立）。
   **むしろこれは現状より良い** — 今は $\Sigma^{\rm MTO}$ の内挿と基底の非周期性が混ざっている。
3. `pwmode = 1` のまま段 2 を回す場合、周期性の破れがどれだけ効くかは測れる（§5 の検証 7）。

## 4. 実装手順

### 4.0 結局どこが置き換わるのか

**`sene` を $H$ に足し込む 1 箇所**、すなわち `getsenex`（`rdsigm2.f90`）の中身だけである。

| | 現状 | MLO 版 |
|---|---|---|
| 内挿 | `call bloch2(qp,ispsigm,sene)` → $\Sigma^{\rm MTO}(q)$（$L\times L$、$L$=230） | $\Sigma^{\rm MLO}(R)$ を同じ対リストで Bloch 和 → $\Sigma^{\rm MLO}(q)$（$M\times M$、$M$=154 or 12） |
| 挟む行列 | `ovliovl` $=\big[(O^{\rm MTO})^{-1}S^{\rm PMT}\big]$ ← 引数 `ovlm` から作る | $\big[(O^{\rm MLO})^{-1}A^\dagger\big]$、$A = S^{\rm PMT}(q)\,z^{\rm MLO}(q)$（式 (15)） |
| 出力 | `senex(ndimh,ndimh)` | 同じ |

出力 `senex(ndimh,ndimh)` は同じなので、呼び出し元の**論理**は変わらない。
ただし現行の引数は `getsenex(qp,isp,ndimh,ovlm)` で **`hamm` を受け取っていない**ので、
$z^{\rm MLO}(k)$ を作るために **`hamm` を 1 つ足す**必要がある:

```
call getsenex(qp, isp, ndimh, ovlm, hamm)   ! hamm = H^LDA(k)、4 箇所とも呼び出し時点で手元にある
```

呼び出しは `m_bandcal.f90`:198,210（通常モード）と `sugw.f90`:429,440（**GW ドライバモード**）の 4 箇所。
どちらのモードも同じ `getsenex` を通るので、**バンドも密度も、GW へ渡す固有関数も、同時に同じ MLO 内挿になる**。
両モードで処方が揃うことは反復の整合に必須である（§3.2.4）。

### 4.1 反復の外（一度だけ）

**$\tilde\chi_\alpha$ は反復ごとに作り直さない。**種 $\chi^{\rm MTO}_{ix(\alpha)}$ も窓のパラメータも固定なのだから、
$\tilde\chi$ を一度決めて全反復で使い回せばよい。これで `mlo` を反復の中で 2 回呼ぶ（サンドイッチする）必要が無くなる。

| # | 実行 | 何に依存するか |
|---|---|---|
| 0a | `lmf --jobgw=0` | **構造と基底だけ**。`m_hamindex0_init()` で `HAMindex0` を書いて即終了する（`main_lmf.f90`:95）。密度にも $\Sigma$ にも依らない |
| 0b | `qg4gw --job=1` | 同上。$q$ メッシュと cutoff から `QGpsi`（offset-Γ 込み）を作るだけ |
| 0c | `lmf --writeham --mlo` | `HamiltonianPMTInfo`（MTO チャネル表の材料） |
| 0d | `mlo --mlo` | **`__mloindex`**: `ndimMTO, ldim, mlomethod, nskip_global` / `ix(1:ndimMTO)` / `fff1, eferm, ecbot`。あわせて `HamRsMLO` |

$z^{\rm MLO}_0(q)$ を**ファイルに置かない**のが実装上の判断である。置くと
「$z^{\rm MLO}_0$ は `__HamiltonianGW` から作る／`__HamiltonianGW` は `lmf --jobgw=1` が書く／
a' はその `lmf --jobgw=1` の中」で循環する。代わりに**索引だけ凍結**し、
$z^{\rm MLO}_0(q)$ は a' の中で $H^{\rm LDA}(q)$ からその場で作り直す（§4.2 a'）。
索引が固定なら処方は同一なので、$\tilde\chi$ は毎回同じ関数になる。

**したがって反復は a3 以降だけでよい。**
現行の `gwsc` は 0a・0b を毎反復走らせているが、どちらも安価なので害は無い（残してよい）。
$\tilde\chi_\alpha$ を 0d で固定するのが要点で、これにより反復の中で MLO を作り直す必要が無くなる。

### 4.2 反復ごと

| 段 | 実行 | 何をする / 出力 |
|---|---|---|
| (a1) | `lmf --jobgw=0` | §4.1 の 0a と同じ。**構造のみ依存なので本来は反復不要** |
| (a2) | `qg4gw --job=1` | §4.1 の 0b と同じ。**同上** |
| **a3** | **`lmf --jobgw=1`**（= `m_sugw_init`） | GW の全 q で固有値問題を解き、`__VxcEvec`（$z^{\psi}$, $V_{xc}$）、`geig`/`cphi`、`__HamiltonianGW` を書く |
| **a'** | **a3 と同じルーチンの中** | $c'(q) = (z^{\psi})^\dagger S^{\rm PMT}(q)\, z^{\rm MLO}_0(q)$ → `__cmlo.data` |
| b | `heftet` / `hbasfp0` / `hvccfp0` / `hsfp0_sc` / `hgw` | $\Sigma^{\psi}$（`SEX2U`, `SEC2U`, …） |
| c | `hqpe_sc` | `sigm`（式 (3)）**と並列に** **`SigmMLO.q`**（式 (11)）。下の注を見よ |
| d | `mlo --mlo` | **`SigRsMLO`**（対称化 → 全 BZ → FFT） |
| e | `lmf` | SCF。`getsenex` が `SigRsMLO` を使う（段 4） |

#### 注: `sigm` と `SigmMLO.q` は並列であって、連鎖ではない

**MTO 基底へ落としてから MLO へ移す、ということはしていない。** `hqpe_sc` の中で、
同一の `se`（= 固有関数基底の $\Sigma^{\psi}-V_{xc}$）から 2 本が**独立に**出る:

```
se (band basis)
 ├─ evec_invt · se · evec_inv  + (窓外)·eseavrmean  →  sigm        230 = MTO   式 (3)
 └─ cmlo^dag  · se · cmlo      + (窓外)·eseavrmean  →  SigmMLO.q   154 = MLO   式 (11)
```

`ev_se_ev`（MTO 枝の中間量）は MLO 枝に一切入らない。
MTO 枝を残してあるのは、現状 `getsenex` が `sigm` を読むからという**過渡的な理由だけ**である。
段 4 が入れば $H$ に実際に入るのは MLO 版になり、`sigm` は従来互換のために書くだけになる。

#### a' を `lmf --jobgw=1` に置く理由

その場に**必要なものが全部揃っている**から:

| 要るもの | a3 の中での正体 |
|---|---|
| $z^{\psi}(q)$ = 現反復の固有ベクトル | `evec`。**`__VxcEvec` に書かれるのと同一の配列**なので、$\Sigma^{\psi}$ が後で使う基底と必ず一致する |
| $S^{\rm PMT}(q)$ | `ovlm`。`hambl` が同じ $q$ で返している |
| $z^{\rm MLO}_0(q)$ | §4.1 で凍結したもの。$H^{\rm LDA}$ で一度だけ作る |
| GW の q リスト（offset-Γ 込み） | `qplist(1:nqirr)`。a3 はこれを回っている |

**基底の一致が「同じ配列を使う」ことで保証される**のが要点である。
`mlo` を別プロセスで走らせて `__HamiltonianGW` から作り直すと、
どの $H$ を書いたかに依存して基底がずれる（§4.1 末尾の実測）。

#### a' の実装  【2026-09-24 実装済み — `sugw.f90`】

`sugw` の $q$ ループの中で、次の順に行う。

1. **`getsenex` の前**に `hamm`（= $H^{\rm LDA}(q)$）と `ovlm`（= $S^{\rm PMT}(q)$）を退避する。
   `zhev_tk4` が両方を破壊するので、退避は必須。
2. 対角化の**後**（`evec` = $z^{\psi}$ が出てから）:
   ```
   call Hreduction(mlomethod, .false., ndimhx, hamm_lda, ovlm_keep, &
                   ndimMTO, ix, fff1, hmo, omo, qp, nev=nxq, zMLO=zm, nskip_auto=nskip)
   sz     = matmul(ovlm_keep, zm)                       ! S^PMT z^MLO_0
   cmlo   = matmul(transpose(dconjg(evec)), sz)         ! c' = z^psi^dag S z^MLO_0
   ```
3. `__cmlo.data` に書く。レコードは `(nbandmx, ndimMTO)`、索引 `isp + nspx*(iq-1)`。
   `__cmlo.info` は最後に master が書く。

`mlomethod`, `ix`, `fff1`, `nskip` は `__mloindex` から読む。
`Hreduction` は窓に `m_readqplist` の `eferm`/`ecbot` を使うが、`sugw` は `qplist.dat` を読まないので
**`set_bandedge(eferm, ecbot)`**（`m_qplist.f90` に新設）で `__mloindex` の凍結値を流し込む。

**`evec` は `__VxcEvec` に書かれるのと同一の配列**である。これにより
$c'$ が $\Sigma^{\psi}$ と同じ band 基底にあることが、規約ではなく**構造的に**保証される。

**制限**: `nspc=2`（SOC）では a' を飛ばす（警告を出す）。スピノル基底での $c'$ は未実装。

これで gwsc から `mlo` の前半呼び出しが消え、**反復ごとの `mlo` は d の 1 回だけ**になった。

### 段 1 — `zMLO` の書き出し  【2026-09-24 完了】

`m_HamPMT.f90` の `HreductionIqibz` ループで `zMLO` を常時要求し、その MTO 行を書き出す。

- `__amlo.data` … 直接アクセス、レコード = `(iq, isp)`、中身 `zMLO(1:ldim, 1:ndimMTO)`（complex(8)）
- `__amlo.info` … `ldim, ndimMTO, nqibz, nspx, recl` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)`

`lmf --writeham --mlo`（= `job_mlo` の第 1 段）で生成される。

### 段 1' — `cmlo` の書き出し  【骨格は実装済み。§4.2 a' の基底そろえが未了 — ここが次の一手】

式 (11) が要るのは `zMLO` ではなく **`cmlo`$_{i\alpha}(q)$**（バンド × MLO、両添字とも $q$ 非依存）。
`Hreduction` の中で作られているが外へ出していない。段 1 と同じ機構で出す。

- `__cmlo.data` … レコード = `(iq, isp)`、中身 `cmlo(1:nband, 1:ndimMTO)`
- `__cmlo.info` … `nband, ndimMTO, nqibz, nspx, recl` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)` /
  **凍結する窓の大域量** `nskip`, $\varepsilon_{\rm cbot}$, $E_F$（任意 $k$ で同じ窓を使うため。§3.2.4）
- **`cmlo` は $H^{\rm LDA}$ から作る**。現状 `HreductionIqibz` は `sigm` があれば $\Sigma$ 込みの $H$ を使うので、
  ここを LDA に固定する（§3.2.4 の循環回避と整合条件）

**注意 — バンド索引の整合**: $\Sigma^{\psi}_{ij}(q)$ は `emax_sigm` 以下のバンド、
`cmlo` は MLO の窓（`nskip` 以上）で張られる。**両者の $i$ の原点と範囲を `info` に明記して突き合わせる**。
窓の外のバンドは式 (11) で重み $\bar\theta$ により**連続に**落ちる（鋭い打ち切りではない）。

### 段 2 — $\Sigma^{\rm MLO}(R)$ の生成  【2026-09-24 実装済み（hqpe_sc + m_HamPMT）】

入力 `sigm`（または $\Sigma^{\psi}$ を直接）+ `__cmlo.data/.info`、出力 `SigRsMLO`。

1. 式 (11) で $\Sigma^{\rm MLO}(q) = (c^{\rm MLO})^\dagger \Sigma^{\psi}(q)\, c^{\rm MLO}$。**MTO を経由しない。**
   $\Sigma^{\psi}$ は `hqpe_sc` が MTO へ落とす前の段にあるので、そこに分岐を入れるのが素直。
2. 既約 $q$ → 全 BZ。**回転則に注意**: $\tilde\chi_\alpha$ は原子中心の $lm$ を種に持つので
   MTO と同じ回転則で回るが、`ix(α)` の対応付けを間違えると全部壊れる。段 3 の検証で必ず捕まる。
3. FFT → $\Sigma^{\rm MLO}(R)$。**格納形式は `HamRsMLO` の `hammr`/`ovlmr` と同一**
   （`npair(ib1,ib2)` / `nlat` / `nqwgt` の対リスト、`(npairmx, ndimMTO, ndimMTO, nspx)`）。
   `HamPMTtoHamRsMLO` の FFT 部分をそのまま流用できる。

### 段 3 — 試験台: `calc_ham_eigen` に足してバンドだけ見る  【2026-09-24 実装済み】

`m_mlo_ham.f90` の `read_ham_rs` で `SigRsMLO` も読み、`calc_ham_eigen` の
`FourierTransform` ブロックで `hammr` と同じ位相で和して `hamm` に加算するだけ。
$O^{\rm MLO}(k)$ は既に同じ場所で作られている。$\Sigma$ を PMT に戻す段は**存在しない**。

これで「$\Sigma$ の内挿を MLO 空間で行ったバンド」が直接得られ、検証 3・4 がすぐ回せる。

### 段 4 — 本体: `getsenex` の差し替え

`rdsigm2.f90` に分岐を入れる。既定は従来どおり。

```
[gw] sigma_mlo = true     # ctrlg のキー（仮）
```

有効時、§4.0 の表のとおり `bloch2` と `ovliovl` を差し替える。任意 $k$ で要るのは
式 (15) の $A(k) = S^{\rm PMT}(k)\,z^{\rm MLO}(k)$ と $O^{\rm MLO}(k)$ の 2 つだけ。

**$z^{\rm MLO}(k)$ の作り直し手順**（内挿しない。その $k$ でゼロから作る）:

1. `__cmlo.info` から **種のリスト `ix(1:ndimMTO)`、`mlomethod`、凍結した窓**
   （`nskip`, $\varepsilon_{\rm cbot}$, $E_F$）を読む。**その場で決め直さない。**
   `Hreduction` は `m_readqplist` の `eferm`, `ecbot` を見るので、そこへ凍結値を入れる。
2. その $k$ で `Hreduction(mlomethod, iprx, ndimh, hamm, ovlm, ndimMTO, ix, fff1, hammout, ovlmout, qp, cmlo, nev, zMLO, nskip)` を呼ぶ。
   中でやっていること:
   1. $H^{\rm LDA}(k)\,c=\varepsilon\,S^{\rm PMT}(k)\,c$ を解く → $\psi^{\rm PMT}_i(k)$, $\varepsilon_i(k)$（`evecpmt`）
   2. MTO ブロックだけを対角化 → $\psi^{\rm MTO}_j(k)$, $\varepsilon^{\rm MTO}_j(k)$（`evecmto`）
   3. `fac`$_{ij}=\langle\psi^{\rm PMT}_i\mid\psi^{\rm MTO}_j\rangle = \big((z^{\psi})^\dagger S^{\rm PMT}_{[:,\,ix]} z^{\rm MTO}\big)_{ij}$
      に窓の重み $\bar\theta_{ij}$ を掛けて `Amat`（$i\le$ `nskip` はゼロ）
   4. `cmlo` $=$ `Amat` $\cdot$ `evecmto`$^\dagger\,S^{\rm MTO}[ix,ix]$
   5. `zMLO` $=$ `evecpmt(:,1:nx)` $\cdot$ `cmlo`
   まとめると、$z^{\psi}$ = `evecpmt`、$z^{\rm MTO}$ = `evecmto` として

$$z^{\rm MLO}(k) = z^{\psi}\,
\Big[\;\bar\theta \odot \big( (z^{\psi})^\dagger\, S^{\rm PMT}_{[:,\,ix]}\, z^{\rm MTO}\big)\Big]\,
(z^{\rm MTO})^\dagger\, S^{\rm PMT}_{[ix,\,ix]}\tag{17}
$$

   （$\odot$ は要素ごとの積）。**材料は 2 つだけ** — **エネルギー因子 $\bar\theta$** と
   **$\langle\chi^{\rm PMT}\mid\chi^{\rm MTO}\rangle = S^{\rm PMT}_{[:,\,ix]}$**（コードの `ovlmx(:,ix)`）。
   固有ベクトルは $H^{\rm LDA}(k)$ を PMT 空間と MTO ブロックで解いて得る。
   **検算**: $\bar\theta\equiv1$ なら完全性から $z^{\rm MLO}$ は種そのもの（$\tilde\chi_\alpha=\chi^{\rm MTO}_{ix(\alpha)}$）に戻る
   — `mlomethod=9` がこの極限である。窓だけが MLO を MTO から隔てている。

3. $A(k) = $ `ovlm` $\cdot$ `zMLO`、$O^{\rm MLO}(k) = $ `zMLO`$^\dagger A$（式 (15)）。
4. 式 (14) で `senex` を組む。

**offset-Γ を含む任意 $k$ で通る**（$\mathrm{FFT}[z]$ を使わないので APW の事情に触れない）。
**循環しない**: MLO は $H^{\rm LDA}$ から作るので、`getsenex` 突入時点で入力が揃っている。
メッシュ点の `cmlo` も同じく LDA 由来にすること（段 1'）。
**コスト**: 2-1 と 2-2 で $k$ あたり対角化 2 回分。

- `mlo_nkabc = n1n2n3` を必須にする。

## 5. 検証手順

| # | 何を | 合格の目安 | 必要な段 |
|---|---|---|---|
| 1 | **$\Sigma^{\rm MLO}(R)$ の減衰** … $\Vert\Sigma(R)\Vert$ vs $\Vert R\Vert$ をチャネル別に。09-24 の MTO 版（Ti3d オンサイト 98.5 → 最近接 6.67 → 以後 0.3〜0.6 meV で頭打ち）と比較 | BvK セル端で頭打ちしない | 段 2 |
| | **【2026-09-25 09:45 暫定合格】** LiTi₂O₄ 6³・76 軌道（LDA 由来、iter 1）。**絶対値では比較できない**: `sigm` は両側に $O^{-1}$ を畳み込んだ量でオンサイト 73.8 eV、$\Sigma^{\rm MLO}=c^\dagger\Sigma^\psi c$ はオンサイト 187 meV。オンサイトで規格化すると $\vert R\vert $=1.4/1.9/2.5/3.5 alat で **MLO 1.8/0.83/0.63/0.36 ×10⁻²** 対 **MTO 4.8/2.9/2.3/2.7 ×10⁻²**。**MTO の方が減衰せず**（5.2 alat まで 2〜3 % で平ら）、MLO は中距離で 2〜3 倍よく落ちる。ただし MLO の対リストは $\vert R\vert \lesssim3.7$ alat で打ち切り。**この測定は $\tilde\chi$ 凍結の入っていないバイナリで取ったので、要再測定**（§6 リスク 3）| | |
| 2 | **メッシュ点での厳密性** … $q$ メッシュ点で $\Sigma^{\rm MLO}$ から作ったバンドを元と比較 | 1 meV 以内 | 段 3 |
| 3 | **補間の滑らかさ** … Γ–X 211 点の残差 rms | 従来 9.8 meV → 6 meV 以下（t2g 模型の実測 3.5 meV が下限の目安） | 段 3 |
| 4 | **模型サイズ** … `mlo_lm` を t2g(12) / 全 EH lm(154) / +EH2 部分追加 と振る | 窓内を 0.5 meV で再現する最小の模型 | 段 3 |
| 5 | **反復の安定性** … 6³ で 10 反復 | 収束解が従来と数 meV 以内、占有 t2g のさざ波（9³ で b34 が 3.1 → 21.9 meV に成長した現象）が成長しない | 段 4 |
| 6 | **他の系での回帰** … Si, NiO, Fe（TestInstall） | 既定（`sigma_mlo` 無効）で完全一致 | 段 4 |
| 7 | **基底の $q$ 周期性の破れ** … $E_n(q)$ と $E_n(q+G_i)$、および $\Sigma^{\rm MLO}(q)$ と $\Sigma^{\rm MLO}(q+G_i)$ を `pwmode` = 1 と 11 で比較（§3.4） | `pwmode=11` で厳密一致。`pwmode=1` での差が内挿の荒れ（〜10 meV）より十分小さいなら 1 のままでよい | 段 2 |

---

## 6. リスク・未解決

1. **バンド索引の突き合わせ** — 段 1'。$\Sigma^{\psi}$ の範囲（`emax_sigm`）と `cmlo` の窓（`nskip` 以上）が
   ずれると式 (11) が静かに壊れる。検証 2（メッシュ点での厳密性）で必ず捕まるので、そこで担保する。
2. **既約 $q$ → 全 BZ の回転** — 段 2-2。`ix(α)` と回転則の対応を誤ると全部壊れる。同じく検証 2 で捕まる。
3. **【重大・2026-09-25 09:50 判明／運用ミス】$\tilde\chi$ 凍結の入っていないバイナリで
   LiTi₂O₄ を全部回していた** — kt1 の `m_sigmlo.f90` に `ZmloRef` が 1 行も無く、
   どの run も `ZmloRef` を書いていなかった。凍結コミット `01b9ea4df` は 09-25 05:25、
   kt1 のソース／バイナリは 05:13。凍結が無いと `getsenex` は呼ばれるたびに
   **その時点の $H$** から $z^{\rm MLO}$ を作り直す一方、$\Sigma^{\rm MLO}(R)$ は
   **前の反復の $\tilde\chi$** で書かれている。式 (14) の $A$ と $\Sigma^{\rm MLO}$ が
   別の基底のものになる。iter 1 は両方 LDA 由来でズレが小さく、**iter 2 が初めての本格的な不整合**
   （`gwdrv` が iter 1 だけ 0 なのと一致）。
   **したがって LiTi₂O₄ の MLO 結果（弦ずれ 23.5→146.5 meV、「Σ が 1.9 倍強い」、iter 10 の局在）は
   すべて保留**。kt1 を同期・再ビルドして LDA から 2 反復で再測定中。

   ここから先は、kt1 で run を投げる前に `strings <bin> | grep -c <新しい文字列>` で
   狙った修正が入っているか確認する。
4. **【解決済み】$\tilde\chi$ が密度とともに動く** — §4.1 は「$\tilde\chi$ はチェーンを通じて
   同じ関数」と書いているが、**実装はそうなっていない**。`getsenex` も a' も $z^{\rm MLO}(k)$ を
   **その時点の $H^{\rm LDA}(k)$** から作り直し、$H^{\rm LDA}$ は密度とともに変わる。一方
   $\Sigma^{\rm MLO}_{\alpha\beta}(R)$ は**前の反復の基底で作られた行列**なので、
   **$\Sigma$ が動く物差しで測り直され続ける**。
   NiO 実測（LDA 密度 → QSGW 1 反復後）: $z^{\rm MLO}$ の相対変化 1.4〜1.9 %、
   軌道ごとの重なり $\langle\tilde\chi^{\rm LDA}\mid\tilde\chi^{\rm QSGW}\rangle$ 最小 0.977〜0.992。
   LiTi₂O₄ 6³ で弦からのずれが 10 反復で 23.5 → 146.5 meV に育つ症状を、これで説明できる。
   `--mlofreeze` は窓（`ix`, `nskip`, `eferm`, `ecbot`）しか凍結しておらず不十分。
   → **`ZmloRef` に $z^{\rm MLO}$ を永続化して凍結した**（commit `01b9ea4df`）。
   なお LiTi₂O₄ の連鎖では $\tilde\chi$ は後半でほぼ収束しており（反復間 0.15〜0.5 %、重なり 0.99998 以上）、
   **荒れが育つ原因ではなかった**。凍結は設計どおりなので実装は残す。
5. **(B) のコスト** — 段 4。`getsenex` の中で $k$ ごとに `Hreduction` を呼ぶので対角化 2 回分増える。
   循環はしない（$H^{\rm LDA}$ から作るので `getsenex` 突入時点で入力が揃っている）。
   重すぎるなら (A) のバンド経路だけを採り、SCF は従来の MTO 内挿のまま、という使い分けもあり得る。
6. **窓の外** — §3.3。E_F+4 eV 以上のバンドは変わる。DOS や光学応答を取るときは注意。
7. **模型の選択** — `mlo_lm` に何を入れるかが事実上のパラメータになる。
   **鉄則（2026-09-24 に NiO で確立）**: **(サイト, lm) ごとに MLO は 1 本まで**。
   MLO は固定した種を band 多様体へ射影したもの（式 (10)）だから、
   窓の中にその性格のバンドが 1 本しか無ければ、種を 2 本（EH と EH2）与えても
   **射影で同じ関数に潰れる**。$O^{\rm MLO}$ が特異に近づき、式 (14) は
   $\big(O^{\rm MLO}\big)^{-1}$ を 2 回挟むので悪化が 2 乗で効く。
   NiO 実測: Ni d に EH2 を足すと $O^{\rm MLO}$ の条件数 $1.3\times10^2\to7.4\times10^4$ に悪化。
   2 本要るなら窓を広げて第 2 のバンドを窓内に入れるしかない（未検証）。
   LiTi₂O₄ では全 EH lm（154）が窓内を 0.5 meV で再現し、荒れは 6.3 meV。
   t2g だけ（12）は荒れ 3.5 meV とより滑らかだが窓が狭い。**既定は全 EH lm** を推奨。
8. **EH2 / PZ** — `mlo_lm2` / `mlo_lm3` で第 2 動径関数・局所軌道も足せるようにしたが、
   **全部（230 チャネル）入れると PMT の完全性を 1 % 以上失って破綻する**（`Hreduction` が停止）。
   部分追加（Ti d だけ 174、O p だけ 178、両方 198）は動く。

---

## 7. 付録 — 今回の実験で分かった「やってはいけないこと」

| やったこと | 結果 |
|---|---|
| Σ(R) の EH2 チャネルの**オフサイトだけ**ゼロ（`ECALJ_SIG_1RAD=2`） | **発散**（E_F = Inf、電子数 638）。近似特異な方向でオンサイトとの打ち消しが壊れる |
| Σ から EH2 ブロックを**丸ごと**射影で落とす（`ECALJ_SIG_EHONLY=2`、Ti 3d のみ） | 計算は通り、**荒れは 9.2 → 5.2 meV**。しかし Σ の値が **平均 292 meV（max 635）動く**。使えない |
| 同上、EH2 全部（`=-1`） | 発散 |
| `nkabc` を GW メッシュに一致 | scf からは補間が消えるが、**バンドは変わらない**（rms 1.5 meV） |
| 空格子球 | QSGW では product basis の作り直しが必要。格子間バンドは `emax_sigm` の窓の外なので割に合わない |

いずれも「**生の MTO 表現では、値を保ったまま自由度を減らせない**」ことを示している。
MLO 表現だけがそれを可能にする、というのが本設計の根拠である。

---

## 8. 実装ログ — 2026-09-24 に踏んだバグと落とし穴

段 1'〜3 を実装して Si と LiTi₂O₄ で回した際に踏んだもの。同じ穴を二度掘らないための記録。

### 8.1 コードのバグ

| # | 症状 | 原因 | 対処 |
|---|---|---|---|
| 1 | $\Sigma$ が一部の $k$ で消え、バンドが X 点で LDA に戻る | `sigmloi` を全 BZ ループの前に **MPI 集約していなかった**。`hammi`/`ovlmi` は `mpibc2_complex` されているのに見落とした | `m_HamPMT` に `mpibc2_complex(sigmloi,...)` を追加（`2ac1552c0`） |
| 2 | `hqpe_sc` が `rotmatMTO` で segfault | `m_mlo_wfs::read_cmlo` は対称操作で回転するが、その `rotmatMTO` は `m_hamindex` の状態を要求する。`hqpe_sc` はそれを用意していない | `get_cmlo_qirr` を新設。`qplistgw` に**そのまま在る** $q$ だけを扱い回転しない（`7039157fa`） |
| 3 | （潜在）`__cmlo.data` のレコードは `(nbandmx, nmlo)` なのに `(nband, nmlo)` の配列へ直接 read していた | `nband < nbandmx` のとき列がずれる。今回は `nband = nbandmx = 322` で顕在化しなかった | 正しい大きさの緩衝を介して読む `read_cmlo_rec` を新設 |
| 4 | メッシュ点で 44〜218 meV ずれる | `__HamiltonianGW` に $H^{\rm LDA}$ を書いたため `cmlo` が $\langle\psi^{\rm LDA}\mid\tilde\chi\rangle$ になり、$\Sigma^{\psi}$（QSGW 固有基底）と**異なる基底を縮約**していた | いったん差し戻し、最終的に a' を `sugw` に実装して解決 |
| 5 | 同 312 meV | 逆の不整合。$\tilde\chi$ を QSGW で作り、$H^{\rm MLO}$ を LDA で作った | 同上 |
| 6 | NiO で EH2 を MLO に足したら伝導帯が 1 eV 動いた | バグではなく**設計の誤り**。同じ (サイト,lm) に 2 本の種を与えると射影で潰れ、$O^{\rm MLO}$ が特異に近づく | §6 リスク 5 の鉄則。EH 1 枚／(サイト,lm) にする |
| 7 | 従来版と MLO 版を**同時に**走らせて比較したら、平均 135 meV・ギャップ 0.87 eV 過小という値が出た | 同時実行による汚染。逐次・クリーンでやり直すと 25 meV・45 meV になった | **比較は必ず逐次・空ディレクトリで。LDA 行が両者で完全一致することを健全性チェックにする** |

### 8.2 運用・手順の落とし穴

| # | 症状 | 原因 |
|---|---|---|
| 6 | バンドが全部 LDA になっていたのに気付かなかった | テストスクリプトが `sigm.<sname>` を消したまま復元せず、lmf が `bndfp (warning): no sigm file found ... LDA calculation only` を出していた。**警告は出るが計算は止まらない** |
| 7 | `mlo` が `gwsc` の中で動かない | `qplist.dat` を必須にしていた（`job_band` 先行が前提）。→ 任意にし、無ければ `eferm` を `efermi.lmf` から取りバンド部を飛ばす |
| 8 | 同上、`HamiltonianPMTInfo` が無い | `lmf --writeham` が先に要る。→ `gwsc` の反復の**外**で 1 回走らせる |
| 9 | `HamRsMLO` に $\Sigma$ を足すと二重計上 | `__HamiltonianPMT` は `getsenex` の**後**に書かれるので既に $\Sigma$ 込み。比較は「LDA の `HamRsMLO` + `SigRsMLO`」対「通常の QSGW バンド」で行う |
| 10 | kt1 のビルドが止まる | `m_mlo_wfs.f90` が nvfortran -O2 で ICE（既知）。該当オブジェクトだけ -O1 で手動コンパイル。`ecaljF_mp`（非 GPU の MP 版）は -O1/-O0 でも通らず未ビルド（kt1 の実行は GPU 版なので支障なし） |

### 8.3 実装は健全だが結果が悪い、を切り分けた手順（2026-09-25）

LiTi₂O₄ で結果が悪かったとき、実装のどこが悪いかを次の順で潰した。診断は環境変数で入る。

| 環境変数 | 何を見るか | LiTi₂O₄/NiO での結果 |
|---|---|---|
| `ECALJ_ZMLO_DUMP=1` | `sugw` の a' と `getsenex` が作る $z^{\rm MLO}$ を両方書き出して比較 | **32 q 点すべてビット単位で一致**（健全） |
| `ECALJ_SIGMLO_RT=1` | メッシュ点で $\Sigma^{\rm MLO}(R)$ の Bloch 和が `__SigmMLO.q` を再現するか | **1.9〜4.5 meV**（最大要素の 0.3 %。健全） |
| `ECALJ_SIGMLO_CHECK=1` | q ごとに MLO 経路と従来経路の `senex` の差 | 全 q で 13〜23 %、offset-Γ が最悪。**q 固有の破綻ではない** |
| （`__cmlo.data` から直接） | $O^{\rm MLO}$ の固有値・条件数 | 36（76 軌道）/ 111（154 軌道）。**健全** |

4 つとも健全だったので、残ったのが §6 リスク 3（$\tilde\chi$ が動く）である。
**「実装が正しくても定式化の前提が破れている」ことがあり得る**という教訓。

### 8.4 測り方の落とし穴

**局所 2 次フィットからの残差は粗さの指標にならない。** バンド交差の位置で必ず尖るので、
交差の多い $t_{2g}$ 多様体では交差だけを拾ってしまう（09-24 の最初の測定はこれで誤判定しかけた）。

**使うべきは「隣接メッシュ点を結ぶ弦からのずれ」**である。$\Sigma$ の内挿誤差はメッシュ点で 0、
その間で最大になるので、この量が内挿のオーバーシュートを直接測る。
バンド分散そのものも含むが、2 つの表現を**同じバンド・同じ区間**で比べる限り差は内挿に由来する。

---

## 9. 基底のセットアップと $\Sigma$ の保持形

2026-09-25 の確認。ここから §13 までが **2 スロット方式**で、§3 の「設計」と §4 の「実装手順」を
いくつかの点で置き換えている。対応は下表。

| §3/§4 の記述 | 置き換えるもの | なぜ |
|---|---|---|
| §10.2 a' が $H^{\rm LDA}$ から $z^{\rm MLO}$ を作る | **§10.2**: 反復 2 以降は $\Sigma$ 入りの $H$ から作る | 標準の MLO（`__HamiltonianPMT`）も senex 加算後の $H$ から作られている |
| `getsenex` が必要に応じて $z$ を作り直す | **§11.1 I1**: 作り直しを禁止し、保存した $z$ でしか読まない | $\Sigma^{\rm MLO}$ の数値は、それを書いた $\tilde\chi$ の刻印を持つ |
| $\tilde\chi$ の置き場は 1 つ | **§10.3**: `ZmloSig` / `ZmloNew` の 2 スロット | ①の中で「読む用 $z_{n-1}$」と「書く用 $z_n$」が同時に要る |
| — | **§11.5**: `nkabc == n1n2n3 == mlo_nkabc` を `gwsc` が事前に要求 | $K_{\rm SCF}\subset K_{\rm GW}$ を成り立たせる |

### 9.1 PMT 基底のセットアップ

そのとおりです。MT 内の**球対称**ポテンシャルから $\phi_l(r,E_\nu)$, $\dot\phi_l$ を作り、
それで smooth Hankel 包絡（MTO）を augment し、APW を足す。
基底は球対称 DFT ポテンシャルだけで決まります。

ただし 2 点:

- **密度が変われば基底も変わる**（弱いが変わる）。
- `pwmode = 11` では APW の集合が $q$ に依存する（$|q+G|$ カット）。

### 9.2 Σ 行列の保持形 — MLO と従来で形が違う

#### MLO 経路

$$\Sigma^{\rm MLO}_{\alpha\beta} = c^\dagger \Sigma^\psi c,
\qquad c = \langle\psi|\tilde\chi\rangle$$

なので、**文字どおり $\langle\tilde\chi_\alpha|\hat\Sigma|\tilde\chi_\beta\rangle$**（行列要素）です。
戻すときは $A = S^{\rm PMT} z^{\rm MLO} = \langle\chi|\tilde\chi\rangle$ を使って

$$\mathrm{senex} = A\,O^{-1}\,\Sigma^{\rm MLO}\,O^{-1}A^\dagger
= \langle\chi|\hat P\hat\Sigma\hat P|\chi\rangle$$

（式 (8)(14)）。整合しています。

#### 従来 `sigm`

MTO ブロックのみ（`ndimsig = nlmto`、APW は入らない）で、しかも**両側に $O^{-1}$ を畳み込んだ形**。
だから今朝オンサイトが 73.8 eV もありました（MLO 側は 187 meV）。
**同じ「⟨基底|Σ|基底⟩」ではありません。**

### 9.3 ⟨ψ|Σ|ψ⟩ から基底表現へ

そのとおりです。ただし窓 $\bar\theta$ の外の状態は落とされ、残りは `eseavrmean` の定数で埋めるので、
**変換が厳密なのは窓の中だけ**です。

---

---

## 10. いつ MLO を作るのか — 1 反復のタイムライン


### 10.0 用語

| | 中身 | 決まる時期 |
|---|---|---|
| **MLO 索引** | `ix`（どの MTO チャネルを種にするか）、`ndimMTO`、窓のパラメータ `eferm`/`ecbot`/`nskip` | **連鎖の最初に 1 度だけ**（`HamRsMLO`）。$\Sigma^{\rm MLO}$ はこのチャネル上の行列なので、`ndimMTO` を途中で変えられない |
| **MLO 係数** $z^{\rm MLO}(q)$ | $\tilde\chi_\alpha(q)=\sum_m\chi^{\rm PMT}_m(q)z_{m\alpha}$ の展開係数 | **毎反復作り直す**（下記 ①） |

以降 $z_n$ = 反復 $n$ で作った係数、$\Sigma_n$ = 反復 $n$ で作った $\Sigma^{\rm MLO}$。

### 10.1 連鎖の最初に 1 度だけ

| 段 | 実行 | 出力 |
|---|---|---|
| 0a | `lmf --jobgw=0` | `HAMindex0`（構造のみ） |
| 0b | `qg4gw --job=1` | `QGpsi`（offset-Γ 込みの q リスト） |
| 0c | `lmf --writeham --mlo` | `HamiltonianPMTInfo` |
| 0d | `mlo --mlo` | **`HamRsMLO`** — 末尾レコードに MLO 索引（`ix`, `ndimMTO`, `mlomethod`, `nskip`, `eferm`, `ecbot`） |

**ここで決まるのは索引だけ**で、$z^{\rm MLO}$ は作りません。

### 10.2 反復 $n$ の中の順序 ← ここが肝心

`gwsc 1` の 1 回がこの順に走ります。

| # | 実行 | $\tilde\chi$ に関して何が起きるか |
|---|---|---|
| ① | **`lmf --jobgw=1`**（[sugw.f90](../../SRC/subroutines/sugw.f90)） | |
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

#### ⑦の注意 — 2 つとも今日ハマった

**(1) `--mlo` を必ず渡す。** `job_band` は `parse_known_args()` で残りの引数をそのまま `lmf` に
渡します（[job_band:91,93](../../SRC/exec/job_band#L91-L93)）。`--mlo` が無いと `c0_mlo` が偽になり、`sigmlo_init` は
冒頭の `if(.not.c0_mlo) return` で即座に戻ります。**MLO 経路は完全に無効**になり、`getsenex` は
従来の `sigm` を使う。`hqpe_sc` は `sigm` と `__SigmMLO.q` を並行して書くので `sigm` は常に
存在し、**黙って従来内挿のバンドが描かれます**。今日の図はすべてこれでした。

**(2) 連鎖ディレクトリの中で走らせない。** `job_band` の第 1 段は `lmf --quit=band` で、
これは **`efermi.lmf` を書き換えます**（[m_bndfp.f90:274](../../SRC/subroutines/m_bndfp.f90#L274)）。MLO の窓は `efermi.lmf` から
`ecbot` オフセットを読むので、連鎖ディレクトリで描くと窓の基準が動きます。
`snap/iter<N>/` を掘って状態一式を置き、その中で走らせる。**そのままスナップショットになります。**

**答え**: MLO（係数 $z$）は **反復ごとに 1 回、①の段 a' で、その反復のハミルトニアンから**作られます。
そして**⑤で初めて「現行」になります**。

#### どの $H$ から作るか — $H^{\rm LDA}$ は初段だけ

$\tilde\chi$ は、これから $\Sigma$ を表現する相手のハミルトニアンの多様体を張っていなければ意味が
ありません。反復 2 以降の $H$ は既に $\Sigma$ を含んでいるので、**$z_n$ も $\Sigma$ 入りの $H$ から作ります**。

| 反復 | a' が使う $H$ |
|---|---|
| 1 | $H^{\rm LDA}$ — まだ $\Sigma$ が無いので（`sigmamode` が偽） |
| 2 以降 | $H = H^{\rm LDA} + \mathrm{senex}(\Sigma_{n-1})$ |

コード上は `zhev_tk4` が `hamm` を破壊するため、senex を足した直後の $H$ を `hamm_qsgw` に
退避して a' で使います（[sugw.f90](../../SRC/subroutines/sugw.f90)）。2 スロット化しているので、`getsenex` 側は `ZmloSig` を
読むだけであり、a' がどの $H$ を使っても読み戻しは壊れません。

**これは標準の MLO 手法と同じです。** [m_bandcal.f90](../../SRC/subroutines/m_bandcal.f90) では

| 行 | 処理 |
|---|---|
| [207](../../SRC/subroutines/m_bandcal.f90#L207) | `hambl` → $H^{\rm LDA}$ |
| [210–211](../../SRC/subroutines/m_bandcal.f90#L210-L211) | `sigmamode` なら `getsenex` → **`hamm = hamm + senex`** |
| [221〜](../../SRC/subroutines/m_bandcal.f90#L221) | `if(writeham)` → **`__HamiltonianPMT` に書く** |

の順で、**`__HamiltonianPMT` には senex を足した後の QSGW ハミルトニアンが入ります**。
`mlo` はこれを読んで `HamRsMLO` を作るので、Si の QSGW バンドなどで使っている標準の MLO は
**元から QSGW の $H$ で作られています**。段 a' が `hamm_lda` を使っていたのが標準から外れていた
方であり、今回の変更はそこを揃えたことになります。

#### なぜ①の `getsenex` が前反復の $z_{n-1}$ を使うのか

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

### 10.3 なぜ 2 スロット要るのか

①の中で $z_{n-1}$（読む用）と $z_n$（書く用）が**同時に生きている**必要があるからです。
1 つしか持てないと、どちらかが必ず食い違います。

| ファイル | 中身 |
|---|---|
| `ZmloSig.<procid>` | 現行 `SigRsMLO` が書かれた $\tilde\chi$。[m_sigmlo.f90 `zmlo_sig_load`](../../SRC/subroutines/m_sigmlo.f90) が読む。追記しない |
| `ZmloNew.<procid>` | 段 a' がこの反復の $H$ から作った $\tilde\chi$（[m_sigmlo.f90 `zmlo_new_append`](../../SRC/subroutines/m_sigmlo.f90)）。⑤で [gwsc](../../SRC/exec/gwsc) が rename |

ランクごとのファイルなのは、段 a' が MPI で $q$ を分担するのと、
`pwmode=11` で `ndimh` が $q$ ごとに変わり固定レコード長が使えないためです。

### 10.4 コード上どこで MLO を作るか

**作る実体は 1 箇所**、[`Hreduction`](../../SRC/subroutines/m_hreduction.f90#L45)（`m_hreduction.f90`）です。
`zMLO` を返します。これを 3 種類の場所から呼んでいます。

| 呼ぶ場所 | いつ | 何を決めるか |
|---|---|---|
| [m_HamPMT.f90](../../SRC/subroutines/m_HamPMT.f90) | 段 0d（`mlo --mlo`）、および MLO バンド | **索引**: `ix`（どの MTO チャネルを種にするか）、`ndimMTO`、窓 `nskip`/`eferm`/`ecbot`。`HamRsMLO` の末尾レコードに保存 |
| [sugw.f90](../../SRC/subroutines/sugw.f90) 段 a'（`NewChiForThisIteration`）| 反復ごと、`lmf --jobgw=1` の中、各 $q$ で | **書く側の $z_n(q)$**。反復 2 以降は `hamm_qsgw`（Σ 入り $H$）、初段は `hamm_lda` から。`ZmloNew.<procid>` へ |
| [m_sigmlo.f90](../../SRC/subroutines/m_sigmlo.f90) `getsenex`（`ChiTildeCache`）| `ZmloSig` に無い $k$ のとき | **読む側の $z(k)$**。その場の $H$ から作り、その `lmf` 実行中だけメモリに保持 |

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

### 10.5 MLO 法の要点 — 入力・作り方・内挿の成立条件

#### user のまとめ（そのとおり）

1. 入力として**メッシュ k 点での $H^{\rm PMT}(k)$**（と $S^{\rm PMT}(k)$）が要る。
2. それがあれば $H^{\rm MTO}(k)$ は部分ブロックとして取り出せる。いま考えている MTO チャネル
   （索引 `ix`）だけを使って、メッシュ k 点で $|{\rm MLO}\rangle$ を $|{\rm PMT}\rangle$ から作る行列
   $z^{\rm MLO}(k)$ が得られる。
3. $\langle {\rm PMT}|X|{\rm PMT}\rangle$ があれば $\langle {\rm MLO}|X|{\rm MLO}\rangle$ に展開でき
   （部分空間にはなる）、これなら FFT で任意の $k$ に内挿できる。$H$ も $S$ も $\Sigma$ も。

コードと突き合わせると 1・2 は完全に一致する。3 も原理としては正しいが、**成立条件が一つある**。

#### 位相（ゲージ）の心配は不要

[m_hreduction.f90:322](../../SRC/subroutines/m_hreduction.f90#L322):

```fortran
cmlo_loc = matmul(Amat, matmul(evecmto^dagger, ovlmx(ix,ix)))
```

コメントどおり $|F^{\rm MLO}_k\rangle = \hat P\,|F^{\rm MTO}_k\rangle$ であり、
$\Psi^{\rm MTO}_j$ も $\Psi^{\rm PMT}_i$ も **$|\Psi\rangle\langle\Psi|$ の形でしか現れない**ので、
対角化の任意位相は相殺する。種は**裸の MTO 基底関数** $\chi^{\rm MTO}_k$（k 非依存の実空間関数の
ブロッホ和）。したがって $\tilde\chi$ にゲージ依存性は無い。

#### 成立条件: $\hat P(k)$ の k についての滑らかさ

$$\hat P(k)=\sum_{ij}|\Psi^{\rm PMT}_i(k)\rangle\,\bar\theta_{ij}(k)\,
\langle\Psi^{\rm PMT}_i(k)|\Psi^{\rm MTO}_j(k)\rangle\langle\Psi^{\rm MTO}_j(k)|$$

$\tilde\chi(k)=\hat P(k)\,\chi^{\rm MTO}(k)$ は、「k 非依存の関数のブロッホ和」を
**k ごとに違う射影で修正したもの**である。q 周期性は保たれるので FFT 自体は正当だが、

> **実空間で短距離になるかどうかは、$\hat P(k)$ が k についてどれだけ滑らかかで決まる。**

では $\hat P$ の k 依存性はどこから来るか。[m_hreduction.f90:308](../../SRC/subroutines/m_hreduction.f90#L308):

```fortran
Amat(i,j) = fac(i,j) * max( fermidist((evl(i)-efrz)/ewfrz), fermidist((evl(i)-ecut)/ewuse) )
```

**窓ではない。** 数字で確認した（user 指摘）:

- **第 1 項は恒等的に 0**。`mlomethod=4` では [`efrz = -1d99`](../../SRC/subroutines/m_hreduction.f90#L246)
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

### 10.6 自己無撞着にどう持っていくか — 採る方式

#### 出発点: $\langle\chi^{\rm PMT}|\tilde\chi\rangle$ は内挿できない

$z^{\rm MLO}$ の入った行列は **PMT 基底の行**を持つ。§3.4 のとおり、PMT のうち
q 周期的なのは **MTO 部分だけ**で、**APW 部分は q 周期的でない**
（`pwmode=11` は $|q+G|$ カットなので G の集合自体が q で変わる）。
行が k ごとに別の関数を指すので、**この行列は実空間に落とせず、内挿できない**。

同じ理由で $H^{\rm PMT}(R)$ を保存して任意の $k$ で $H^{\rm PMT}(k)$ を再構成する道も無い。
**$H^{\rm PMT}(k)$ は、その密度で実際にハミルトニアンを組む以外に入手できない。**

→ したがって **$\tilde\chi$ は毎回作るしかない。そして作るにはハミルトニアンが要る。**

#### 採る方式

> **反復 $n$ で $\Sigma_n$ を作るとき、その場の $H$（＝反復 $n$ 開始時点のハミルトニアン）から
> $z^{\rm MLO}(k)$ を「後で読むことになる全ての $k$」で作り、保存する。
> 以後 $\Sigma_n$ を読むときは必ずその $z$ を使う。**

$\Sigma^{\rm MLO}(R)$ 自体は MLO 基底上の行列なので実空間に落とせ、任意の $k$ に内挿できる。
戻すための $A(k)=S^{\rm PMT}(k)\,z^{\rm MLO}(k)$ は**保存済みのもの**を使う。
これで**内挿できないものを内挿しない**という条件を満たす。

#### 必要な k の網羅 — ここが現状との差

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

#### これで何が保証されるか

- $\Sigma_n$ は**常に、それを作った $H$ の $\tilde\chi$** で読み戻される（任意の $k$ で）。
- $\tilde\chi$ は**反復ごとに更新される**（窓が古いバンド位置に張り付かない）。
- 内挿するのは $\Sigma^{\rm MLO}(R)$ だけ。$\langle\chi^{\rm PMT}|\tilde\chi\rangle$ は内挿しない。

#### 残る近似

反復 $n$ の最後の `lmf` は SCF で密度を動かすので、その中で $H$ は変わる。
$z$ は反復開始時の $H$ のものを使い続けるため、SCF ループ内では厳密には整合しない。
ただし
- 同じ $z$ を使い続けるので **$\Sigma$ の読み方は一貫している**（揺れない）
- 自己無撞着に達すれば密度が止まるので、この近似は消える

つまり **途中の反復でのみ残る近似**であり、収束を妨げる種類のものではない、というのが見立て。

#### 実装量

- `sugw` 段 a' を「q リスト」から「必要な k 集合」へ拡張（k の列挙と、そこでの `Hreduction` 呼び出し）
- SCF の k メッシュとバンド経路を `lmf --jobgw=1` の時点で知る必要がある
  （前者は `nkabc` から、後者は `syml.<sname>` から）
- 少し面倒だが、機構自体（`ZmloSig`/`ZmloNew` と昇格）は既にある

---

---

## 11. 新方式 — 実装仕様


§10.6 の結論を実装できる形にしたもの。**これを実装する。**

### 11.1 不変条件

| | 内容 |
|---|---|
| **I1** | $\Sigma^{\rm MLO}_n$ は**必ず $z_n$ で読み戻す**。読む $k$ が何であっても例外を作らない |
| **I2** | $z_n$ は**$\Sigma_n$ を作ったハミルトニアン**から作る（反復 $n$ 開始時点の $H$）|
| **I3** | **内挿するのは $\Sigma^{\rm MLO}(R)$ だけ**。$\langle\chi^{\rm PMT}\vert \tilde\chi\rangle$ は内挿しない（APW 行が q 周期的でないため不可能）|
| **I4** | MLO **索引**（`ix`, `ndimMTO`）は連鎖を通じて固定。$\Sigma^{\rm MLO}$ がそのチャネル上の行列だから |

I1 と I2 を両立させるために、$z$ を**後で読む全ての $k$** で先に作っておく。これが新方式の核心。

### 11.2 k 集合 $K$

$$K = K_{\rm GW}\ \cup\ K_{\rm SCF}\ \cup\ K_{\rm band}$$

| | 中身 | 出どころ |
|---|---|---|
| $K_{\rm GW}$ | GW が要求する q（既約 BZ ＋ offset-Γ）| `QGpsi` / `m_qplist` |
| $K_{\rm SCF}$ | SCF の k メッシュ | `nkabc`（`n1n2n3` と一致するとは限らない）|
| $K_{\rm band}$ | バンド経路上の k | `syml.<sname>` から生成 |

**$K$ は `lmf --jobgw=1` の中で列挙できなければならない。**
$K_{\rm band}$ のためにバンド経路を事前に固定する必要があるが、診断用の固定経路なので問題ない。
`syml` が無ければ $K_{\rm band}=\emptyset$ とし、バンドは描かない（描くなら `syml` を置く）。

#### 11.2.1 【後回し】MLO 空間で描けば $K_{\rm band}$ は不要になる

バンドを **PMT 空間**で描く限り、経路上の各 $k$ で $A(k)=S^{\rm PMT}(k)z^{\rm MLO}(k)$ が要るので
$K_{\rm band}$ を先に回しておく必要がある。

しかし **MLO 空間で描けばその必要が無い**。$H^{\rm MLO}(R)$ も $\Sigma^{\rm MLO}(R)$ も
**どちらも MLO 基底上の行列で q 周期的**なので、両方を任意の $k$ にフーリエ内挿して

$$\bigl[H^{\rm MLO}(k)+\Sigma^{\rm MLO}(k)\bigr]\,c = \varepsilon\,O^{\rm MLO}(k)\,c$$

を $n_{\rm MLO}\times n_{\rm MLO}$ で解けばよい。**PMT への引き上げが一切要らない。**

- $K = K_{\rm GW}\cup K_{\rm SCF}$ だけになる（$K_{\rm band}$ 不要）
- §11.5 の 3 点目（`syml` を `--jobgw=1` の前に用意）も不要になる
- 機構は既にある（`job_mlo`、`HamRsMLO` の `hammr`/`ovlmr`）

描かれるのは「縮約模型のバンド」であって PMT ハミルトニアンのバンドそのものではないが、
**MLO 模型が $\Sigma$ を再現できているかを見るのが目的**なので、診断としてはむしろ素直。

**2026-09-25 時点では後回し**（user 判断）。まず $K=K_{\rm GW}\cup K_{\rm SCF}$ で実装し、
バンドは従来どおり PMT 空間で描く（そのぶん経路上は作り直しになる）。

#### 11.2.2 PMT 空間のバンドはどうするか — 3 段階

「フルの PMT でバンドを描くのを諦めるのか」への答えは、段階によって違う。

| | 状況 | 基底の整合 | 方針 |
|---|---|---|---|
| **収束後** | 密度が止まっている | あとから描いても $z$ は $z_n$ と一致する | **問題なし。発表用のバンドは損なわれない** |
| **反復の途中** | 密度が $\Sigma$ を作った時点から動いている | メッシュ点は `ZmloSig`、その間は作り直し ＝ **混在** | **診断には使わない**。代わりに MLO 空間で描く（§11.2.1）|
| **途中でも PMT が欲しいとき** | 同上 | — | **可能**。下記 |

**3 段目は実装できる。** `sugw` は各 q で自分で `hambl` を呼んでいる
（[sugw.f90:492-508](../../SRC/subroutines/sugw.f90#L492-L508)）。つまり GW ドライバの中で
バンド経路の $k$ についても `hambl` → `getsenex` → `Hreduction` → `zmlo_new_append` を
回せば、**$\Sigma$ を作ったのと同じ $H$ で**経路上の $z$ が得られる。
機構は同じ subroutine 内に揃っているので「a' の第二ループを足す」程度の作業
（k リストの MPI 配布と `ndimh` の可変性には注意）。

**当面の方針**: 2 段目を採る（MLO 空間で描く）。必要になれば 3 段目を足す。

### 11.3 ファイル

| ファイル | 中身 | 寿命 |
|---|---|---|
| `SigRsMLO` | $\Sigma^{\rm MLO}(R)$ ＋ 対リスト ＋ `ix`, `ib_tableM` | 反復をまたぐ |
| `ZmloSig.<procid>` | **$k\in K$ 全部**での $z_n(k)$。現行 `SigRsMLO` と対 | 反復をまたぐ |
| `ZmloNew.<procid>` | 作成中の $z_{n+1}(k)$ | 反復内 |
| `HamRsMLO` 末尾 | `ix`, `ndimMTO`, `mlomethod`, `nskip`, `eferm`, `ecbot` | 連鎖を通じて固定 |

`ZmloSig` の被覆が $K_{\rm GW}$ だけから $K$ 全体に広がるのが、現行との唯一の構造差。

### 11.4 1 反復の手順

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

#### 11.4.1 ①の `getsenex` の式

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

これが $H(k)$ に足される。式 (12)(14)(15) に対応。

#### なぜ読む側は $z$ を作れないのか（整合以前の問題）

`getsenex(qp, isp, ndimh, ovlm, hamm)` に渡る `hamm` は **senex を足す前の $H$**（$\Sigma$ 抜き）である。
`getsenex` はまさにその senex を**これから計算する**ために呼ばれているので、
**$\Sigma$ 入りの $H$ はその時点で存在しない**。循環している。

対して `sugw` の段 a' は**対角化の後**に走るので `hamm_qsgw`（$=H+\mathrm{senex}$）が既にある。

| | $\Sigma$ 入りの $H$ | $z$ をどうするか |
|---|---|---|
| **書く側**（段 a'）| **ある**（senex 加算後）| その場で作れる |
| **読む側**（`getsenex`）| **無い**（これから作る）| **作れない → 保存したものを読むしかない** |

したがって 2 スロットは「整合を保つための工夫」ではなく、**読む側に選択肢が無いことの帰結**である。
旧実装が `getsenex` で $z$ を作り直していたのは $\Sigma$ 抜きの $H$ で代用していたということで、
**定義からして別の $\tilde\chi$** を作っていたことになる。

**新方式で変わるのはここではない。** 式は同じで、**$z^{\rm MLO}(k)$ をどこから持ってくるか**だけが変わる
（その場で作り直す → $\Sigma$ を作った $H$ で作って保存したものを引く）。

### 11.5 現行コードとの差分と実装状況（2026-09-25 時点）

| | 内容 | 状態 |
|---|---|---|
| 1 | **`m_sigmlo`** — $k$ が `ZmloSig` に無かったら**黙って作り直さず停止**（`ECALJ_MLO_ALLOW_REBUILD=1` で従来動作）| **実装済み** `fb9662fa0` |
| 2 | **`gwsc`** — `nkabc == n1n2n3 == mlo_nkabc` を**計算前に**確認し、違えば直し方つきで停止 | **実装済み** `fb9662fa0` |
| 3 | $K_{\rm band}$ の被覆 | **未実装**。バンドは依然 PMT 空間で作り直し。§11.2.1（MLO 空間で描く）で解消する方針 |
| 4 | 昇格を `gwsc` から `mlo` 側へ移す | **未実装**。`SigRsMLO` 書き出しと不可分にすると取りこぼしが構造的に消える |

**$K_{\rm SCF}$ の拡張は不要だった。** 当初は a' のループを $K$ 全体へ広げるつもりだったが、
任意の $k$ で $z$ を作るにはその $k$ で $H(k)$ を組む必要があり（`hambl` の複製）侵襲的。
代わりに **`nkabc == n1n2n3` を要求**すれば $K_{\rm SCF}\subset K_{\rm GW}$ が成り立つ。
QSGW では密度に $\Sigma$ の内挿を入れないため元々メッシュを揃えて使うので、実害のない制約。

なお Fortran 側では $K_{\rm SCF}$ を確認できない（GW ドライバでは `bz_nkp = 0`、
SCF の k リストが構築されていない）。そのため `gwsc` で `ctrlg` を読んで確認している。

#### GaAs での実測（`Samples/GetStarted/GaAs`、2 原子、`pwmode=11`）

| 設定 | 結果 |
|---|---|
| `nkabc`=4³, `n1n2n3`=2³（不一致）| SCF が $q=(0,0,0.5)$ で外れる。`ZmloSig` の 7 点に無い。**停止する** |
| `nkabc`=`n1n2n3`=2³（一致）| **2 反復完走**。GW ドライバ `loaded=7 / MISS=0`、最後の `lmf` `loaded=3 / MISS=0`、昇格 2 回 |

`ZmloSig` の中身は 2³ の既約 BZ（Γ, L, X）＋ offset-Γ 4 点 = 7 レコード。
外れた $q=(0,0,0.5)$ は 4³ 固有の点で、`ndimh` の不一致ではなく**本当に q が無い**ケースだった。

**テストに使えなかった系**: NiO は AF（`AF=1/-1`）で SCF が重く、さらに `hqpe_sc` の
`if(laf) exit` でスピン 2 の $\Sigma^{\rm MLO}$ が 0 になる懸念がある（**未確認の問題**）。
`Samples/MLOsamples/C` は MLO 専用サンプルで `[product_basis]` が GW 用に整っておらず
`hbasfp0 --job=3` が segfault する。

### 11.6 検証（これが揃って初めて「新方式で測った」と言える）

| 見るもの | 期待 |
|---|---|
| `llmfgw01` / `llmf` / `llmf_band` の「作り直し」回数 | **0**（1 回でも出たら $K$ の列挙漏れ）|
| `gwsc: promoted N ZmloNew.* -> ZmloSig.*` | 毎反復 |
| `loaded ZmloSig ... records=` | 各実行で $\vert K\vert $ 相当 |
| `ECALJ_ZMLO_DUMP` | a' と `getsenex` の $z$ が一致 |
| `ECALJ_SIGMLO_RT` | $\Sigma^{\rm MLO}(q)\to R\to q$ の往復誤差 |

### 11.7 残る近似

反復 $n$ の最後の `lmf` は SCF で密度を動かすが、$z$ は反復開始時のものを使い続ける。
したがって SCF ループ内では $z$ と $H$ が厳密には対応しない。ただし

- 同じ $z$ を使い続けるので **$\Sigma$ の読み方は一貫している**（k でも SCF ステップでも揺れない）
- 自己無撞着に達すれば密度が止まるので、この近似は消える

**途中の反復にのみ残る近似**であり、収束を妨げる種類のものではない、というのが見立て。

### 11.8 コスト

$K_{\rm band}$ ぶんの `Hreduction` が増える。1 k あたり PMT の対角化 1 回（$O(n_{\rm dimh}^3)$）。
211 点の経路なら GW ドライバの実行時間に数分乗る程度で、GW 諸段に比べれば小さい。

---

---

## 12. Q&A


### Q1. §11.4.1 の $O^{-1}$ は「部分空間の逆」か「PMT 空間の逆」か。どちらがよいか

**候補は 2 つ。**

| | 形 | 次元 |
|---|---|---|
| **(a) 部分空間の逆**（現実装）| $O^{\rm MLO}=\langle\tilde\chi\vert \tilde\chi\rangle=z^\dagger S z$、$\hat P=\vert \tilde\chi\rangle (O^{\rm MLO})^{-1}\langle\tilde\chi\vert $ | $n_{\rm MLO}\times n_{\rm MLO}$ |
| **(b) PMT 空間の逆** | $S^{-1}$ を使う形 | $n_{\rm dimh}\times n_{\rm dimh}$ |

**MLO については (a) しかない。** 理由は 2 つ。

1. $\tilde\chi_\alpha$ は **PMT 基底関数そのものではなく線型結合**である。
   従来の MTO 経路では部分空間が PMT の**部分ブロック**なので「$S$ の小行列を逆にする」ことに
   意味があり、実際 [getsenex](../../SRC/subroutines/rdsigm2.f90#L67) は `ovlm(1:ndimsig,1:ndimsig)` を
   逆にしている。**MLO には逆にすべき小行列が存在しない。**
2. $S$ を計量とする非直交基底で $\hat P^2=\hat P$ を満たす射影子は
   $|\tilde\chi\rangle O^{-1}\langle\tilde\chi|$ に**一意に決まる**。
   $S^{-1}$ を使う形は、$\Sigma^{\rm MLO}$（$n_{\rm MLO}^2$ 個の拘束）から
   PMT 係数 $\sigma$（$n_{\rm dimh}^2$ 個）を決めることになり**劣決定**で、
   最小性条件を課すと結局 (a) に戻る。

**未確認の点（正直に）**: 従来経路で `hqpe_sc` が書く `sene` の規約。
[main_hqpe.sc.f90](../../SRC/subroutines/main_hqpe.sc.f90) は**擬似逆**
$\psi^+=(\psi^\dagger\psi)^{-1}\psi^\dagger$ で挟み、`getsenex` 側は **$S_{\rm sub}^{-1}$** で挟み返す。
この 2 つが整合して $\langle\chi|\hat\Sigma|\chi\rangle$ を再現するのか、代数的に追い切れていない
（オンサイトが 73.8 eV という大きさは、行列要素そのものではなく双対側の量であることを示唆する）。

**数値で決着できる**: `ECALJ_SIGMLO_CHECK=1` は各 q で MLO 版と MTO 版の `senex` を直接比較する。
規約が食い違っていれば**系統的な因子**として現れる。
2026-09-25 に走らせたときは「全 q で一様に 13〜23 %」だったが、
あれは基底不整合のあるバイナリでの測定なので、**新方式で測り直す**（§11.6 の検証項目に追加）。

### Q2. なぜ従来法にはこの面倒が無いのか

引き上げ行列 $A$ が**その場の重なり行列の切り出し**だから。
従来の部分空間は PMT 基底の部分ブロックそのものなので $A=S^{\rm (full,sub)}$ であり、
保存すべきものが何も無い。MLO は $A=S^{\rm PMT}(k)\,z^{\rm MLO}(k)$ で、
$z$ が **$k$ 依存かつ $H$ 依存**である。ここが構造的な非対称。

### Q3. $z^{\rm MLO}$ を実空間に落として内挿すればよいのでは

**できない。** $z^{\rm MLO}$ の入った行列は **PMT 基底の行**を持ち、
PMT のうち q 周期的なのは MTO 部分だけで **APW 部分は q 周期的でない**
（`pwmode=11` は $|q+G|$ カットなので G の集合自体が q で変わる）。
行が $k$ ごとに別の関数を指すので実空間に落とせない。
同じ理由で $H^{\rm PMT}(R)$ を保存して $H^{\rm PMT}(k)$ を再構成する道も無く、
**$H^{\rm PMT}(k)$ はその密度で実際にハミルトニアンを組む以外に入手できない**。
→ だから §11 の方式（読む $k$ を全部先に回しておく）になる。

### Q4. 1 反復で MLO を何回作っているのか

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

### Q5. 窓 $\bar\theta$ を広げれば $\Sigma(R)$ は短距離になるか

**ならない。** 2026-09-25 に「窓が狭いのでは」と書いたが撤回した。

- `mlomethod=4` では第 1 項が恒等的に 0（[`efrz = -1d99`](../../SRC/subroutines/m_hreduction.f90#L246)、
  「active only for mlomethod=3」）。効くのは第 2 項だけ。
- 第 2 項の位置は $\varepsilon^{\rm cut}_j=\max(\varepsilon_{\rm cbot}+\Delta,\ \varepsilon^{\rm MTO}_j)$。
  MTO のみの固有値は PMT のそれより上にあるので通常は $\varepsilon^{\rm MTO}_j$ が選ばれ、
  引数は $(\varepsilon^{\rm PMT}_i-\varepsilon^{\rm MTO}_j)/w$ という**バンド差**になって $k$ について滑らか。
- 幅も狭くない。LiTi₂O₄ は `mlo_w = 2.0` eV、`mlo_delta = 2.0` eV。

$\hat P(k)$ の $k$ 依存性は窓ではなく、
**固有ベクトル $|\Psi^{\rm PMT}_i(k)\rangle\langle\Psi^{\rm PMT}_i(k)|$ 自体が $k$ で回ること**から来る。

### Q6. $\tilde\chi$ にゲージ（位相）の任意性は入らないか

**入らない。** [m_hreduction.f90:322](../../SRC/subroutines/m_hreduction.f90#L322) で
$|F^{\rm MLO}_k\rangle=\hat P\,|F^{\rm MTO}_k\rangle$ の形になっており、
$\Psi^{\rm MTO}_j$ も $\Psi^{\rm PMT}_i$ も $|\Psi\rangle\langle\Psi|$ の形でしか現れないので
対角化の任意位相は相殺する。種は**裸の MTO 基底関数** $\chi^{\rm MTO}_k$（$k$ 非依存の実空間関数の
ブロッホ和）。

---

---

## 13. TODO


優先順に。番号は依存関係ではなく、いまの見立ての順序。

### A. すぐやるもの

| | 内容 | なぜ | 状態 |
|---|---|---|---|
| **A1** | **kt1 を同期・再ビルド**して LiTi₂O₄ 6³ を新方式で 3 反復 | `nkabc = n1n2n3 = mlo_nkabc = 6³` なので新チェックは通る。今日の比較対象（従来 194→86→32 meV、旧 MLO 415→158→193 meV）と同じ長さで並べられる | kt1 は `d0ca0aba8`（保護なし）。同期 1 分＋ビルド 6 分＋LDA 3 分＋16 分×3 |
| **A2** | **バンドを MLO 空間で描く**（§11.2.1）| $H^{\rm MLO}(R)$ も $\Sigma^{\rm MLO}(R)$ も q 周期的なので任意の $k$ で厳密。今日ずっと判定を濁らせた**基底混在が構造的に消える** | **GaAs では動作、LiTi₂O₄ で異常終了**（下記 A2'）|
| **A2'** | 上の異常終了を直す | `job_mlo` の第 1 段は `lmf --writeham --mkprocar --noinv --mlo` で、**`--noinv`** のため q 集合が段 0c のものと違う。`--mlofreeze` で保持している `HamRsMLO` の対リスト（`npair`/`nlat`、段 0c の q 集合で作成）と食い違い、`HamPMTtoHamRsMLO` の `iqiloop 3/16` で `double free or corruption (out)`。np=1 でも同じなので MPI ではない。GaAs（q が少ない）では通った | **未着手** |
| **A3** | **昇格を `gwsc` から `mlo` 側へ移す** | `SigRsMLO` の書き出しと不可分になり、今日の「`tbin_frozen/gwsc` が古くて一度も走らなかった」型の取りこぼしが構造的に消える | 数十行 |

**A1 の判定基準**: $E_{\rm HF}$ と反復間の変化で見る。バンドの荒れは A2 が入るまで参考扱い。

### B. 未確認・要検証

| | 内容 | なぜ気になるか |
|---|---|---|
| **B1** | **AF（`laf`）で $\Sigma^{\rm MLO}$ のスピン 2 が 0 になっていないか** | `main_hqpe.sc.f90` の `if(laf) exit` でスピンループを抜ける。MLO 分岐に AF の扱いが見当たらない。NiO 系が全滅する可能性 |
| **B2** | **`mlo --mlofreeze --mlo` が LDA の $H$ から MLO を作り直している件**（§12 Q4）| `__HamiltonianPMT` は段 0c の LDA のまま。$\Sigma$ の対称化・回転は索引レベルなので汚していないはずだが**未確認** |
| **B3** | **従来経路の $\Sigma$ の規約整合**（§12 Q1）| `hqpe_sc` は擬似逆、`getsenex` は $S_{\rm sub}^{-1}$。`ECALJ_SIGMLO_CHECK=1` で MLO 版と MTO 版の `senex` を比べれば系統因子として出るはず |
| **B4** | **窓（`eferm`/`ecbot`）が LDA 固定のまま** | `Hreduction` は QSGW 固有値を LDA の $E_F$ と比べている。正しい直し方は判明済み: `sugw` 側だけで `call set_bandedge(eferm, eferm + (ecbot_a - eferm_a))`。一度誤実装で iter 1 を壊して revert 済み |
| **B5** | **初回（LDA から）の混合を半歩にするか否か** | `mixsigma` は履歴が無いと $x_0=0$ から線形混合するので、iteration 1 の `sigm` と $\Sigma^{\rm MLO}$ は $\beta\Sigma^{\rm out}$ になる（2026-09-25 v8 で 812.37 → 406.18 を確認）。user の想定は「初回は混合しない」。変えるなら従来チェーンも取り直し |
| **B6** | **`gwsc N`（N≥2）で `__SigmMLO.q.prev` が退避されない** | MLO の後片付け（`.prev` への退避を含む）が反復ループの外。`gwsc 1` を繰り返す運用では無害だが、`gwsc N` かつ `ECALJ_MLO_MIX=1` では 2 反復目以降 `hqpe_sc` が古い $\Sigma$ を $x_0$ に読む。直し方: 後片付けを関数にしてループ先頭でも呼ぶ（冪等）|

### C. あとで

| | 内容 |
|---|---|
| **C1** | 途中反復でも PMT バンドが欲しくなったら、a' に第二ループを足す（§11.2.2 の 3 段目）|
| **C2** | 154 軌道での再測定（`--mlo` 付き）。部分空間の広さの効果。kr5 に入力一式を用意済み（`-np 8` に修正済み）|
| **C3** | kr5 の `gwsc` を通しで動かす（`-np 8` で再開できるはず。LDA は kt1 と完全一致を確認済み）|
| **C4** | 収束後に PMT バンドを描いて従来法と比較（収束後なら基底は整合する）|

### D. 運用（今日の教訓）

| | 内容 |
|---|---|
| **D1** | ラン投入前に **`sync_ecalj_src.sh --check-all`**（ソース刻印）|
| **D2** | ビルド後に **`strings <lib> \| grep -c <新機能の文字列>`**（実体は `.so` 側。実行ファイルは thin wrapper）|
| **D3** | **`ls -la <tbin>/gwsc`** が symlink か実体コピーか |
| **D4** | ラン開始直後にログの印（`promoted`, `loaded ZmloSig`, `MISS`, `mloON`）|

今日はこの 4 点を怠って 3 回無駄にした（ソース未同期 / `FC` 取り違えでビルド失敗して古い `.so` が残留 / `tbin` が前日の実体コピー）。**3 回とも見た目は正常**だった。
