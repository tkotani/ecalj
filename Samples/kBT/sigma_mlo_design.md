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

### 3.0 記号 — 基底・双対基底・MLO・固有関数を区別する

| 記号 | 意味 | 次元 / 索引 |
|---|---|---|
| $\chi^{\rm PMT}_m$ | **PMT 基底関数**（MTO + APW） | $m = 1\ldots n_{\rm dimh}$、**APW の本数は $k$ に依存** |
| $\chi^{\rm MTO}_\mu$ | その MTO 部分 | $\mu = 1\ldots L$（`ldim` = 230）、**$k$ 非依存** |
| $\chi^{\rm MTO}_{ix(\alpha)}$ | **MLO の種**。`ix(`$\alpha$`)` 番目の MTO 基底関数（固定） | 新しい記号は要らない |
| $\tilde\chi_\alpha = P\,\chi^{\rm MTO}_{ix(\alpha)}$ | **MLO**（コードの `F^MLO`）。band 多様体への射影 | $\alpha = 1\ldots M$（`ndimMTO` = 154） |
| $\psi^{\rm PMT}_i$ | **固有関数**（基底ではない） | |
| `cmlo`$_{ik} = \langle \psi^{\rm PMT}_i\mid\tilde\chi_\alpha\rangle$ | MLO の**固有関数**基底での係数 | |
| `zMLO`$_{m\alpha}$ | MLO の **PMT 基底**での係数、$\tilde\chi_\alpha=\sum_m \chi^{\rm PMT}_m z^{\rm MLO}_{m\alpha}$ | $(n_{\rm dimh}, M)$ |
| $P_a$ | **PAW（augmentation）チャネル** φ, φ̇ | $a = 1\ldots N_a$（`ndima` = 700）、**$k$ 非依存** |
| $B_{am}(k)$ | PMT 基底関数の augmentation 係数 | lmf が任意 $k$ で厳密に作る |

コードの `F^MTO_k` は $\chi^{\rm MTO}_{ix(\alpha)}$、`F^MLO_k` は $\tilde\chi_\alpha$ に対応する。
**MLO の索引には $\alpha,\beta$ を使い、波数 $k,q$ と区別する。**

$\Sigma$ は `sigm` に $\Sigma^{\rm MTO}_{\mu\nu}(q)=\langle\chi^{\rm MTO}_\mu|\hat\Sigma|\chi^{\rm MTO}_\nu\rangle$ として入っている（$L\times L$）。

### 3.1 変えられるのは部分空間だけ

§2.5–2.6 のとおり、展開の枠組み（$\big(O\big)^{-1}$ を両側に挟む非直交の射影）は既に一般的な形であり、
**APW ブロックも埋まっている**。したがって設計の自由度は

> **どの部分空間で $\Sigma$ を保持し、内挿するか**

の一点に尽きる。現状は MTO（$L = 230$、非直交、EH/EH2 が準線形従属、裾が長い）。
以下、MLO（§3.2）と PAW チャネル（§3.3）を検討する。
なお **どちらも MTO を経由しない** — 式 (1) の $\Sigma^{\psi}$ から直接その部分空間の行列を作る。

### 3.2 MLO を部分空間に使う【本命】

#### 3.2.1 MLO ↔ PMT の変換行列は既に得ている（段 1 実装済み）

まず確認しておくべき事実として、**MLO を PMT 基底へ展開する係数行列はコードが既に作っている**。
`Hreduction`（`SRC/subroutines/m_hreduction.f90`）が引数 `zMLO` で返すのがそれで、

$$\tilde\chi_{\alpha q} \;=\; \sum_m \chi^{\rm PMT}_{mq}\; z^{\rm MLO}_{m\alpha}(q),
\qquad z^{\rm MLO} \in \mathbb{C}^{\,n_{\rm dimh}(q)\times M}\tag{9}
$$

二段で組まれている。第一段は、**固定した MTO 種** $\chi^{\rm MTO}_{ix(\alpha)}$ を
band 多様体へ射影した係数（固有関数基底、コードの `cmlo`）

$$c^{\rm MLO}_{i\alpha}(q) \;=\; \sum_j
\big\langle \psi^{\rm PMT}_{iq}\mid\psi^{\rm MTO}_{jq}\big\rangle^{*}\,
\big\langle \psi^{\rm MTO}_{jq}\mid\chi^{\rm MTO}_{ix(\alpha)}\big\rangle\tag{10}
$$

第二段は、固有ベクトル $z^{\psi}_{mi}(q)$（`evecpmt`、$|\psi_{iq}\rangle=\sum_m|\chi^{\rm PMT}_{mq}\rangle z^{\psi}_{mi}$）を掛けて PMT 基底へ戻す:

$$z^{\rm MLO}(q) \;=\; z^{\psi}(q)\, c^{\rm MLO}(q)\tag{11}
$$

この構成の性質（内挿の可否を左右する）:

| 性質 | 内容 |
|---|---|
| **ゲージ** | 種 $\chi^{\rm MTO}_{ix(\alpha)}$ が $q$ に依らず固定なので **射影ゲージ**。固有ベクトルの任意位相・縮退内の任意回転は式 (10) で相殺し、$z^{\rm MLO}(q)$ は $q$ の滑らかな関数になる（Wannier の gauge fixing に相当する作業が不要） |
| **直交性** | 既定では**非直交**。$O^{\rm MLO}(q) = (z^{\rm MLO})^\dagger S^{\rm PMT}(q)\, z^{\rm MLO}$。`--mlo_ortho` のときだけ Löwdin 直交化（`MLOLowdinOrthogonalization`）、`c0_mlo_diagnorm` で対角規格化。いずれも任意 |
| **次元** | $M$ = `ndimMTO`（LiTi2O4 全 MTO で 154、t2g モデルで 12）。行数 $n_{\rm dimh}(q)$ は MTO 230 + APW（$q$ 依存） |
| **出力** | `m_HamPMT.f90` が `__amlo.data` / `__amlo.info` に direct access で書き出す（1 record = 1 個の $(q,\sigma)$）。`info` は `ldim, ndimMTO, nqibz, nspx, mrecamlo` と種の索引 `ix(1:ndimMTO)`、`qibz(1:3,1:nqibz)` |

**ただし現状のファイルは MTO 行だけ**（`writem(..., zMLO(1:ldim,1:ndimMTO))`、$230\times154$）で、
APW 行は捨てている。これは §3.2.3 の障害そのものである。

#### 3.2.2 式 (3)–(5) は MTO を MLO に置き換えるだけで済む

§2 の流れのうち **(3) ψ→基底 / (4) FFT / (5) Bloch 和** は、MTO を $\tilde\chi$ に
読み替えるだけで成立する。しかも式 (3) の MLO 版は **`zMLO` を経由しない**:

$$\Sigma^{\rm MLO}_{\alpha\beta}(q) \;=\; \langle \tilde\chi_\alpha|\hat\Sigma(q)|\tilde\chi_\beta\rangle
\;=\; \sum_{ij} \big(c^{\rm MLO}_{i\alpha}\big)^{*}\,\Sigma^{\psi}_{ij}(q)\,c^{\rm MLO}_{j\beta}
\;=\;\big((c^{\rm MLO})^\dagger\,\Sigma^{\psi}\,c^{\rm MLO}\big)_{kl}\tag{12}
$$

$c^{\rm MLO}_{i\alpha}=\langle\psi^{\rm PMT}_{iq}\mid\tilde\chi_{\alpha q}\rangle$（式 (10)）の第一添字は
**バンド指標**なので、APW の本数が $q$ に依ることは**ここには入らない**。
$\psi$ は正規直交だから逆行列も不要。$\Sigma^{\psi}$ から一発で取れる。

続く 2 段はそのまま:

$$\Sigma^{\rm MLO}_{\alpha\beta}(R) = \frac{1}{N_q}\sum_q \Sigma^{\rm MLO}_{\alpha\beta}(q)e^{-iqR},
\qquad
\Sigma^{\rm MLO}_{\alpha\beta}(k) = \sum_R \Sigma^{\rm MLO}_{\alpha\beta}(R)\,\overline{e^{ikR}}\tag{13}
$$

$\alpha,\beta$ は MLO の通し番号（$M$ = `ndimMTO`）で **$q$ 非依存**だから、FFT も Bloch 和も
MTO 版と同一の実装で回る。**内挿される自由度は $L^2 n_k^3 = 230^2\cdot n_k^3$ から
$M^2 n_k^3 = 154^2\cdot n_k^3$（t2g モデルなら $12^2\cdot n_k^3$）に減り、
しかも準線形従属な方向が取り除かれている** — §1 の表の問題に直接効く。

#### 3.2.3 式 (6)(7): $\tilde\chi$ の実空間表現を経由する

残るのは、任意 $k$ で

$$\big[\hat\Sigma\big]_{mn}(k) = \sum_{klk'l'} A_{m\alpha}(k)\,
\big(O^{\rm MLO}\big)^{-1}_{kk'}\,\Sigma^{\rm MLO}_{k'l'}(k)\,
\big(O^{\rm MLO}\big)^{-1}_{l'l}\,A^{*}_{n\beta}(k),
\qquad A_{m\alpha}(k) = \big\langle \chi^{\rm PMT}_{mk}\mid\tilde\chi_{\alpha k}\big\rangle\tag{14}
$$

を作る段（式 (6)(7) の MLO 版）。必要なのは $A(k)$ と $O^{\rm MLO}(k)=\langle\tilde\chi_\alpha|\tilde\chi_\alpha\rangle$ の 2 つだけ。

**やってはいけない道**: $z^{\rm MLO}_{m\alpha}(q)$ を要素ごとにフーリエ内挿すること。
行索引 $m$ は PMT 基底で、APW の本数が $k$ ごとに違う（$|k+G|<$ cutoff の $G$ 集合が変わる）ので
$C(R)=\mathrm{FFT}[z^{\rm MLO}(q)]$ は定義できない。

**正しい道**: $\tilde\chi$ は**実空間の局在関数**であり、その表現を経由する。
$$\tilde\chi_{k}(\mathbf r - \mathbf R) \;=\; \frac{1}{N_q}\sum_{q} e^{-i q R}\,\tilde\chi_{\alpha q}(\mathbf r),
\qquad
\tilde\chi_{\alpha q}(\mathbf r) = \sum_m \chi^{\rm PMT}_{mq}(\mathbf r)\, z^{\rm MLO}_{m\alpha}(q)\tag{15}
$$
この $\tilde\chi_\alpha(\mathbf r)$ は BvK 超格子上で厳密に決まる **$k$ に依らない 1 個の関数**である。中身は 2 部分:

| 部分 | 表現 | 任意 $k$ での扱い |
|---|---|---|
| MTO 部分 | $\sum_{\mu R'} d_{\mu\alpha}(R')\,\chi^{\rm MTO}_\mu(\mathbf r-\mathbf R')$、$d = \mathrm{FFT}\big[z^{\rm MLO}_{\rm MTO}\big]$ — 索引 $\mu$ は $q$ 非依存なので FFT できる（`__amlo.data` にあるのはこれ） | lmf が既に持つ「Bloch 和した MTO どうしの重なり」の機構でそのまま計算できる |
| 平滑（APW）部分 | 超格子の平面波係数 $z^{\rm MLO}_{(qG)\alpha}$ の集合。局在関数なのでフーリエ成分は連続 | 局在領域に切って $\hat f(\mathbf k+\mathbf G)$ を求め直す（リサンプル） |

そして式 (14) が要求するのは $\tilde\chi$ を **$k$ の基底で展開すること**ではなく、
**$\chi^{\rm PMT}_{mk}$ との重なり積分**である:

$$A_{m\alpha}(k) \;=\; \sum_{R} e^{i k R}\,
\big\langle \chi^{\rm PMT}_{mk}\,\big|\,\tilde\chi_\alpha(\mathbf r-\mathbf R)\big\rangle\tag{16}
$$

これは $\tilde\chi$ が $k$ の基底の張る空間に入っているかどうかとは無関係に、常に定義される。
$O^{\rm MLO}(k)$ も同じ形（$\chi^{\rm PMT}$ を $\tilde\chi$ に替えるだけ）。
$\Sigma^{\rm MLO}(R)$ を実空間で切って式 (13) で Bloch 和するのと、**同じ種類の操作**である。

**残る近似**は、局在関数の裾を BvK セルで切ることだけ。
これは現状の式 (4)(5)（$\Sigma^{\rm MTO}(R)$ を BvK セルで切る）が既に受け入れている近似と同質で、
しかも $\tilde\chi$ は $\Sigma$ より局在が良いはずだから条件は緩い。

**実装上いちばん楽な分解**は、球内を PAW チャネル $P_a$（$k$ 非依存、$B(k)$ は lmf が任意 $k$ で厳密に作る）、
球外の平滑部分を超格子平面波、とすること。すなわち §3.3 の道具立てを**戻す段にだけ**使う。
§3.2.2（どの部分空間で $\Sigma$ を保持するか = MLO）と §3.3（どう戻すか = PAW）は**併用できる**。

### 3.3 PAW（augmentation）チャネル【部分空間にも、戻す道具にも使える】

**出発点は §2.1 の式 (1)(2) そのもの**。`hqpe_sc` がそれを式 (3) で MTO 基底に落とすのが現状だが、
式 (3) は**途中産物**にすぎない。$\hat\Sigma$ を別の部分空間で保持したいなら、
**MTO を経由せず $\Sigma^{\psi}$ から直接**取ればよい。

$P_a$（MT 球内の φ, φ̇）は **$k$ 非依存の固定索引**であり、$B(k)$ は **lmf が任意 $k$ で厳密に作る**。
したがって $\Sigma$ をこのチャネルで保持すれば、任意 $k$ への戻しが逆行列なしで厳密に書ける。
また §3.2.3 のとおり、**$\Sigma$ の部分空間は MLO のまま、球内の表現だけ $P_a$ にする**という併用もできる。

**(1) GW の出力から直接**

固有関数の augmentation 係数 `cphi` は $\psi_i = \sum_a P_a\,c_{ai}$（球内）なので
$\langle P_a|\psi_i\rangle = (\Pi c)_{ai}$（$\Pi$ = 球内の重なり `ppj`）。演算子 $\hat\Sigma = \sum_{ij}|\psi_i\rangle\Sigma^\psi_{ij}\langle\psi_j|$ の PAW 行列は

$$\boxed{\ \Sigma^{\rm PAW}(q) \;=\; (\Pi\,c)\;\Sigma^{\psi}(q)\;(\Pi\,c)^\dagger\ }\tag{17}
$$

`hqpe_sc` が今 式 (3) を作っているところを、これに差し替える。
`cphi` も `ppj` も GW 側に既にある。

**(2) 内挿**

$\Sigma^{\rm PAW}(q)$ を FFT して $\Sigma^{\rm PAW}(R)$、任意 $k$ で Bloch 和。**内挿されるのはここだけ。**

**(3) 戻す**

同じ形で書くと、$\Pi_{ab} = \langle P_a\mid P_b\rangle$（`ppj`）を使って

$$\big[\hat\Sigma\big]_{mn} = \sum_{acdb} \langle\chi^{\rm PMT}_m|P_a\rangle\,
\big(\Pi\big)^{-1}_{ac}\,\Sigma^{\rm PAW}_{cd}(k)\,\big(\Pi\big)^{-1}_{db}\,
\langle P_b|\chi^{\rm PMT}_n\rangle\tag{18}
$$

ここで $\langle\chi^{\rm PMT}_m|P_a\rangle = (B^\dagger \Pi)_{ma}$ なので $\Pi^{-1}$ が両側で約分して

$$\boxed{\ \Sigma^{\rm PMT}_{mn}(k) \;=\; \big(B^\dagger(k)\,\Sigma^{\rm PAW}(k)\,B(k)\big)_{mn}\ }\tag{19}
$$

**$B(k)$ は厳密**、逆行列も不要。`getsenex` は 1 行になる。

**局在性**: 部分波は球内に厳密に閉じているので、MTO 包絡（smooth Hankel、裾が長い）より
$\Sigma(R)$ の減衰が速いはず。09-24 の測定で Ti 3d ブロックが BvK セル端で頭打ちだったのは、
包絡の裾が主因の可能性が高い。**要検証**（§5 の検証 1）。

**次元**: `ndima` = 700（`lmxa` = 4）。l を切れば

| 切り方 | チャネル数 |
|---|---|
| 全部（lmxa = 4） | 700 |
| l ≤ 3（= `lmx`） | 448 |
| l ≤ 2 | 252 |
| 原子別（Ti は d まで、O は p まで、Li は s,p） | 〜150 |

**どこで切るかは実測で決める**（§5 の検証 2）。

### 3.4 何を失うか

PAW チャネルの外 = **MT 球の外（格子間）の $\Sigma$**。
現状の MTO 表現も球外を MTO 包絡の裾で表しているだけなので、優劣は自明ではない。**検証 3 で測る。**

## 4. 実装手順

### 段 1 — $A(q)$ の書き出し  【2026-09-24 完了 — 実装済みはここまで】

`m_HamPMT.f90` の `HreductionIqibz` ループで `zMLO` を常時要求し、その MTO 行を書き出す。

- `__amlo.data` … 直接アクセス、レコード = `(iq, isp)`、中身 `zMLO(1:ldim, 1:ndimMTO)`（complex(8)）
- `__amlo.info` … `ldim, ndimMTO, nqibz, nspx, recl` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)`

既存の SOC 用 `zMLO` をそのまま使うので追加コストはほぼゼロ。
`lmf --writeham --mlo`（= `job_mlo` の第 1 段）で生成される。

### 段 2 — $\Sigma^{\rm PAW}(R)$ の生成

新規ツール（または `hqpe_sc` の後段）。入力 `sigm`, `__amlo.data/.info`、出力 `SigRsMLO`。

1. `hqpe_sc` で $\Sigma^{\rm PAW}(q) = (\Pi c)\,\Sigma^{\psi}(q)\,(\Pi c)^\dagger$ を作り、
   `sigm` の代わりに（または併記して）書き出す。`cphi`（= $c$）も `ppj`（= $\Pi$）も GW 側に既にある。
   **MTO 基底は経由しない。**
2. 既約 q → 全 BZ へ展開。**$P_a$ は原子中心の球面調和なので回転則が明快**
   （MTO/MLO のように「行と列で回転則が違う」問題が無い）。
3. FFT → $\Sigma^{\rm PAW}(R)$、WS 最短ベクトルの対リストで保存。

### 段 3 — `getsenex` の差し替え

`rdsigm2.f90` に分岐を入れる。既定は従来どおり。

```
[gw] sigma_mlo = true     # ctrlg のキー（仮）
```

有効時は §3.5 の式で $\Sigma^{\rm PMT}(k)$ を作る。必要な量:
- $\Sigma^{\rm PAW}(k)$ … $\Sigma^{\rm PAW}(R)$ の Bloch 和
- $B(k)$ … augmentation 係数。**lmf が任意 $k$ で厳密に作る**（内挿不要、循環無し）

### 段 4 — gwsc への組み込み

PAW チャネル版は $B(k)$ が厳密なので **反復の中で特別な手当ては要らない**。
MLO を重ねる場合（§3.6）だけ、直前の反復の `zMLO` を使う（`mlo_nkabc = n1n2n3` を必須）。

---

## 5. 検証手順

| # | 何を | 合格の目安 |
|---|---|---|
| 1 | **$\Sigma^{\rm PAW}(R)$ の減衰** … $\Vert\Sigma(R)\Vert$ vs $\Vert R\Vert$ をチャネル別に | MTO 表現（09-24 の測定）より速く減衰。BvK セル端で頭打ちしない |
| 2 | **次元の切り詰め** … l の上限を 4/3/2 と変えて $\Sigma$ の値と荒れを見る | 値が動かない範囲で最小の l |
| 3 | **メッシュ点での厳密性** … §3.3(3) の $\Sigma^{\rm PMT}(k)$ を q メッシュ点で元と比較 | バンドで 1 meV 以内 |
| 4 | **補間の滑らかさ** … Γ–X 211 点の残差 rms | 従来 9.8 meV → 6 meV 以下 |
| 5 | **反復の安定性** … 6³ で 10 反復 | 収束解が従来と数 meV 以内、占有 t2g のさざ波が成長しない |
| 6 | **他の系での回帰** … Si, NiO, Fe（TestInstall） | 既定（`sigma_mlo` 無効）で完全一致 |

---

## 6. リスク・未解決

1. **ゲージ** — §3.5(b)。検証 1 が最優先。駄目なら段 3 は $A(k)$ を LDA から作る案に落ちるが、
   メッシュ点での厳密性（検証 3）を失う。
2. **既約 q → 全 BZ の回転** — 段 2 の注意点。$A$ の回転則を誤ると全部壊れる。
   検証 3 で必ず捕まるので、そこで担保する。
3. **窓の外** — §3.6。E_F+4 eV 以上のバンドは変わる。DOS や光学応答を取るときは注意。
4. **模型の選択** — `mlo_lm` に何を入れるかが事実上のパラメータになる。
   LiTi₂O₄ では全 EH lm（154）が窓内を 0.5 meV で再現し、荒れは 6.3 meV。
   t2g だけ（12）は荒れ 3.5 meV とより滑らかだが窓が狭い。**既定は全 EH lm** を推奨。
5. **EH2 / PZ** — `mlo_lm2` / `mlo_lm3` で第 2 動径関数・局所軌道も足せるようにしたが、
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
