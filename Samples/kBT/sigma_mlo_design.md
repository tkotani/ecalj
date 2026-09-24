# 自己エネルギーの MLO 表現による内挿 — 設計書

2026-09-24 起草。対象: ecalj 開発者。
背景データは [kBT_research.md](kBT_research.md) の 2026-09-23〜24 のエントリ（図表番号はそちらを参照）。

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
3. **(B) のコスト** — 段 4。`getsenex` の中で $k$ ごとに `Hreduction` を呼ぶので対角化 2 回分増える。
   循環はしない（$H^{\rm LDA}$ から作るので `getsenex` 突入時点で入力が揃っている）。
   重すぎるなら (A) のバンド経路だけを採り、SCF は従来の MTO 内挿のまま、という使い分けもあり得る。
4. **窓の外** — §3.3。E_F+4 eV 以上のバンドは変わる。DOS や光学応答を取るときは注意。
5. **模型の選択** — `mlo_lm` に何を入れるかが事実上のパラメータになる。
   LiTi₂O₄ では全 EH lm（154）が窓内を 0.5 meV で再現し、荒れは 6.3 meV。
   t2g だけ（12）は荒れ 3.5 meV とより滑らかだが窓が狭い。**既定は全 EH lm** を推奨。
6. **EH2 / PZ** — `mlo_lm2` / `mlo_lm3` で第 2 動径関数・局所軌道も足せるようにしたが、
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

### 8.2 運用・手順の落とし穴

| # | 症状 | 原因 |
|---|---|---|
| 6 | バンドが全部 LDA になっていたのに気付かなかった | テストスクリプトが `sigm.<sname>` を消したまま復元せず、lmf が `bndfp (warning): no sigm file found ... LDA calculation only` を出していた。**警告は出るが計算は止まらない** |
| 7 | `mlo` が `gwsc` の中で動かない | `qplist.dat` を必須にしていた（`job_band` 先行が前提）。→ 任意にし、無ければ `eferm` を `efermi.lmf` から取りバンド部を飛ばす |
| 8 | 同上、`HamiltonianPMTInfo` が無い | `lmf --writeham` が先に要る。→ `gwsc` の反復の**外**で 1 回走らせる |
| 9 | `HamRsMLO` に $\Sigma$ を足すと二重計上 | `__HamiltonianPMT` は `getsenex` の**後**に書かれるので既に $\Sigma$ 込み。比較は「LDA の `HamRsMLO` + `SigRsMLO`」対「通常の QSGW バンド」で行う |
| 10 | kt1 のビルドが止まる | `m_mlo_wfs.f90` が nvfortran -O2 で ICE（既知）。該当オブジェクトだけ -O1 で手動コンパイル。`ecaljF_mp`（非 GPU の MP 版）は -O1/-O0 でも通らず未ビルド（kt1 の実行は GPU 版なので支障なし） |

### 8.3 測り方の落とし穴

**局所 2 次フィットからの残差は粗さの指標にならない。** バンド交差の位置で必ず尖るので、
交差の多い $t_{2g}$ 多様体では交差だけを拾ってしまう（09-24 の最初の測定はこれで誤判定しかけた）。

**使うべきは「隣接メッシュ点を結ぶ弦からのずれ」**である。$\Sigma$ の内挿誤差はメッシュ点で 0、
その間で最大になるので、この量が内挿のオーバーシュートを直接測る。
バンド分散そのものも含むが、2 つの表現を**同じバンド・同じ区間**で比べる限り差は内挿に由来する。
