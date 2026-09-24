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

以下、そこへ至る段を並べる。段 3 は lmf に手を入れずに数字を先に見るための**試験台**で、
本体は段 4 である。

### 段 1 — `zMLO` の書き出し  【2026-09-24 完了】

`m_HamPMT.f90` の `HreductionIqibz` ループで `zMLO` を常時要求し、その MTO 行を書き出す。

- `__amlo.data` … 直接アクセス、レコード = `(iq, isp)`、中身 `zMLO(1:ldim, 1:ndimMTO)`（complex(8)）
- `__amlo.info` … `ldim, ndimMTO, nqibz, nspx, recl` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)`

`lmf --writeham --mlo`（= `job_mlo` の第 1 段）で生成される。

### 段 1' — `cmlo` の書き出し  【未着手・ここが次の一手】

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

### 段 2 — $\Sigma^{\rm MLO}(R)$ の生成

入力 `sigm`（または $\Sigma^{\psi}$ を直接）+ `__cmlo.data/.info`、出力 `SigRsMLO`。

1. 式 (11) で $\Sigma^{\rm MLO}(q) = (c^{\rm MLO})^\dagger \Sigma^{\psi}(q)\, c^{\rm MLO}$。**MTO を経由しない。**
   $\Sigma^{\psi}$ は `hqpe_sc` が MTO へ落とす前の段にあるので、そこに分岐を入れるのが素直。
2. 既約 $q$ → 全 BZ。**回転則に注意**: $\tilde\chi_\alpha$ は原子中心の $lm$ を種に持つので
   MTO と同じ回転則で回るが、`ix(α)` の対応付けを間違えると全部壊れる。段 3 の検証で必ず捕まる。
3. FFT → $\Sigma^{\rm MLO}(R)$。**格納形式は `HamRsMLO` の `hammr`/`ovlmr` と同一**
   （`npair(ib1,ib2)` / `nlat` / `nqwgt` の対リスト、`(npairmx, ndimMTO, ndimMTO, nspx)`）。
   `HamPMTtoHamRsMLO` の FFT 部分をそのまま流用できる。

### 段 3 — 試験台: `calc_ham_eigen` に足してバンドだけ見る  【lmf 無改造。先にこれをやる】

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
