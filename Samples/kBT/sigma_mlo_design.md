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
| $\chi^{\rm PMT}_m$ | **PMT 基底関数**（MTO + APW） | $m = 1\ldots n_{\rm dimh}$、**APW の本数は $k$ に依存** |
| $\chi^{\rm MTO}_\mu$ | その MTO 部分 | $\mu = 1\ldots L$（`ldim` = 230）、**$k$ 非依存** |
| $\chi^{\rm MTO}_{ix(\alpha)}$ | **MLO の種**。`ix(`$\alpha$`)` 番目の MTO 基底関数（固定） | 新しい記号は要らない |
| $\tilde\chi_\alpha$ | **MLO**（コードの `F^MLO`）。種をエネルギー窓付きで band 多様体へ射影したもの（式 (10)(11)）。**冪等な射影ではない** | $\alpha = 1\ldots M$（`ndimMTO` = 154） |
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

#### 3.2.1 MLO ↔ PMT の変換行列は既に得ている（段 1 実装済み）

まず確認しておくべき事実として、**MLO を PMT 基底へ展開する係数行列はコードが既に作っている**。
`Hreduction`（`SRC/subroutines/m_hreduction.f90`）が引数 `zMLO` で返すのがそれで、

$$\tilde\chi_{\alpha q} \;=\; \sum_m \chi^{\rm PMT}_{mq}\; z^{\rm MLO}_{m\alpha}(q),
\qquad z^{\rm MLO} \in \mathbb{C}^{\,n_{\rm dimh}(q)\times M}\tag{9}
$$

二段で組まれている。第一段は、**固定した MTO 種** $\chi^{\rm MTO}_{ix(\alpha)}$ を
band 多様体へ**エネルギー窓付きで**射影した係数（固有関数基底、コードの `cmlo`）:

$$c^{\rm MLO}_{i\alpha}(q) \;=\; \sum_j \bar\theta_{ij}(q)\;
\big\langle \psi^{\rm PMT}_{iq}\mid\psi^{\rm MTO}_{jq}\big\rangle\,
\big\langle \psi^{\rm MTO}_{jq}\mid\chi^{\rm MTO}_{ix(\alpha)}\big\rangle,
\qquad c^{\rm MLO}_{i\alpha}=0 \ \ (i \le \texttt{nskip})\tag{10}
$$

ここで $\psi^{\rm MTO}_{jq}$ は MTO ブロックだけを対角化した固有関数、
$\bar\theta$ が**窓の重み**（コードの `Amat` $=$ `fac` $\times\,\bar\theta$）:

$$\bar\theta_{ij}(q) \;=\; \max\Big[\,
f\!\Big(\frac{\varepsilon_i(q)-\varepsilon_{\rm frz}}{w_{\rm frz}}\Big),\;
f\!\Big(\frac{\varepsilon_i(q)-\varepsilon^{\rm cut}_{j}}{w}\Big)\Big],
\qquad f(x)=\frac{1}{e^{x}+1}\tag{11}
$$

$\varepsilon^{\rm cut}_j = \max(\varepsilon_{\rm cbot}+\Delta,\ \varepsilon^{\rm MTO}_j(q))$（既定の `mlomethod=4`）。
パラメータは `mlo_delta` $=\Delta$、`mlo_wfrz` $=w_{\rm frz}$、`mlo_w` $=w$、`mlo_down`。
床 $\varepsilon_{\rm cbot}+\Delta$ は **$q$ 非依存**に取る（$k$ ごとに動かすと窓の意味が $k$ ごとに変わり、金属で悪化する）。

**これは冪等な射影ではない**。$\bar\theta<1$ の重みが入り（`--mlo_gs` の Gram–Schmidt を使えばなおさら）、
$\tilde\chi_\alpha$ は「窓の中のバンドで書ける範囲で種に最も近い局在関数」である。
**窓の外は鋭く切られるのではなく、幅 $w$ で連続に落ちる。**

第二段は、固有ベクトル $z^{\psi}_{mi}(q)$（`evecpmt`、$|\psi_{iq}\rangle=\sum_m|\chi^{\rm PMT}_{mq}\rangle z^{\psi}_{mi}$）を掛けて PMT 基底へ戻す:

$$z^{\rm MLO}(q) \;=\; z^{\psi}(q)\, c^{\rm MLO}(q)\tag{12}
$$

この構成の性質（内挿の可否を左右する）:

| 性質 | 内容 |
|---|---|
| **窓** | $\bar\theta$（式 (11)）で窓の中のバンドだけが寄与する。$\Sigma^{\rm MLO}$ は**窓の中でのみ $\hat\Sigma$ を忠実に表す** — これは欠陥ではなく設計意図そのもの（§3.3） |
| **ゲージ** | 種 $\chi^{\rm MTO}_{ix(\alpha)}$ が $q$ に依らず固定なので **射影ゲージ**。固有ベクトルの任意位相・縮退内の任意回転は式 (10) で相殺し、$z^{\rm MLO}(q)$ は $q$ の滑らかな関数になる（Wannier の gauge fixing に相当する作業が不要） |
| **直交性** | 既定では**非直交**。$O^{\rm MLO}(q) = (z^{\rm MLO})^\dagger S^{\rm PMT}(q)\, z^{\rm MLO}$。`--mlo_ortho` のときだけ Löwdin 直交化（`MLOLowdinOrthogonalization`）、`c0_mlo_diagnorm` で対角規格化。いずれも任意 |
| **次元** | $M$ = `ndimMTO`（LiTi2O4 全 MTO で 154、t2g モデルで 12）。行数 $n_{\rm dimh}(q)$ は MTO 230 + APW（$q$ 依存） |
| **出力** | `m_HamPMT.f90` が `__amlo.data` / `__amlo.info` に direct access で書き出す（1 record = 1 個の $(q,\sigma)$）。`info` は `ldim, ndimMTO, nqibz, nspx, mrecamlo` と種の索引 `ix(1:ndimMTO)`、`qibz(1:3,1:nqibz)` |

**ただし、この先で実際に要るのは `cmlo` の方である**（式 (13)）。
現状 `__amlo.data` が持っているのは `zMLO` の MTO 行だけ（`writem(..., zMLO(1:ldim,1:ndimMTO))`、$230\times154$）で、
`cmlo` は `Hreduction` の中で作られてそのまま捨てられている。→ 段 1' で書き出す（§4）。

#### 3.2.2 式 (3)–(5) は MTO を MLO に置き換えるだけで済む

§2 の流れのうち **(3) ψ→基底 / (4) FFT / (5) Bloch 和** は、MTO を $\tilde\chi$ に
読み替えるだけで成立する。しかも式 (3) の MLO 版は **`zMLO` を経由しない**:

$$\Sigma^{\rm MLO}_{\alpha\beta}(q) \;=\; \langle \tilde\chi_\alpha|\hat\Sigma(q)|\tilde\chi_\beta\rangle
\;=\; \sum_{ij} \big(c^{\rm MLO}_{i\alpha}\big)^{*}\,\Sigma^{\psi}_{ij}(q)\,c^{\rm MLO}_{j\beta}
\;=\;\big((c^{\rm MLO})^\dagger\,\Sigma^{\psi}\,c^{\rm MLO}\big)_{\alpha\beta}\tag{13}
$$

$c^{\rm MLO}_{i\alpha}=\langle\psi^{\rm PMT}_{iq}\mid\tilde\chi_{\alpha q}\rangle$（式 (10)）の第一添字は
**バンド指標**なので、APW の本数が $q$ に依ることは**ここには入らない**。
$\psi$ は正規直交だから逆行列も不要。$\Sigma^{\psi}$ から一発で取れる。

続く 2 段はそのまま:

$$\Sigma^{\rm MLO}_{\alpha\beta}(R) = \frac{1}{N_q}\sum_q \Sigma^{\rm MLO}_{\alpha\beta}(q)e^{-iqR},
\qquad
\Sigma^{\rm MLO}_{\alpha\beta}(k) = \sum_R \Sigma^{\rm MLO}_{\alpha\beta}(R)\,\overline{e^{ikR}}\tag{14}
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
式 (14) は新規実装ですらなく、既存の経路に 1 本足すだけである。$O^{\rm MLO}(k)$ も既に任意 $k$ で得られている。

**「APW の本数が $k$ 依存」が効くのは 1 箇所だけ**: $z^{\rm MLO}_{m\alpha}(q)$ を**要素ごとにフーリエ内挿**しようとしたときである。
行索引 $m$ が PMT 基底なので $\mathrm{FFT}[z^{\rm MLO}]$ は定義できない。**その手を使わなければよい**だけで、
$\tilde\chi$ 自体は BvK 超格子上の局在関数として厳密に決まっている:

$$\tilde\chi_{\alpha}(\mathbf r-\mathbf R) = \frac{1}{N_q}\sum_{q} e^{-i q R}\sum_m \chi^{\rm PMT}_{mq}(\mathbf r)\, z^{\rm MLO}_{m\alpha}(q)\tag{15}
$$

#### 3.2.4 使い方は 2 通り

**(A) バンド / QP エネルギーが目的なら、追加実装はほぼ無い。**
$H^{\rm MLO}(R)$ に $\Sigma^{\rm MLO}(R)$ を足して `calc_ham_eigen` に渡すだけ。
$\Sigma$ を PMT 基底に戻す必要は無いので、式 (6)(7) に相当する段自体が存在しない。
09-24 の測定で、MLO 表現の Γ–X 粗さは MTO 内挿の 9.8 meV に対し 3.5 meV（t2g 12 軌道）/ 6.3 meV（154 軌道）、
メッシュ点の値は 0.1–0.5 meV で再現していた。

**(B) SCF（電子密度）まで回すなら、$\Sigma$ を PMT ハミルトニアンに戻す必要がある。**
lmf は PMT 基底で密度を作るので、任意 $k$ で

$$A_{m\alpha}(k) = \big\langle\chi^{\rm PMT}_{mk}\mid\tilde\chi_{\alpha k}\big\rangle = \big(S^{\rm PMT}(k)\,z^{\rm MLO}(k)\big)_{m\alpha}\tag{16}
$$

が要る。**これも内挿しない**: $z^{\rm MLO}(k)$ は式 (10)–(12) の射影をその $k$ で実行すれば直接得られる
（`Hreduction` は `qplist.dat` の任意 $k$ リストで走る。`job_mlo` がバンドプロットで実際にやっていること）。
コストは $k$ ごとの対角化 1 回分。$\mathrm{FFT}[z]$ を使わないので APW の本数が $k$ 依存でも構わない。

**整合条件**: $\Sigma^{\rm MLO}(q)$ を作ったときの band 多様体と、$A(k)$ を作るときの多様体が
**同じ処方**でなければならない（例: どちらも「前反復の $H$ の固有関数」）。
種 $\chi^{\rm MTO}_{ix(\alpha)}$ が固定なのでゲージは自動的に揃う（§3.2.1 の表）。

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

インタフェースは変わらない。呼び出し元（`m_bandcal.f90`、`sugw.f90` の 2 箇所）は無改造で、
**バンドも密度も SCF も同時に MLO 内挿になる**。

新規に要る入力は $z^{\rm MLO}(q)$ ただ 1 つ。$S^{\rm PMT}(q)$ は既に引数 `ovlm` で渡っている。

以下、そこへ至る段を並べる。段 3 は lmf に手を入れずに数字を先に見るための**試験台**で、
本体は段 4 である。

### 段 1 — `zMLO` の書き出し  【2026-09-24 完了】

`m_HamPMT.f90` の `HreductionIqibz` ループで `zMLO` を常時要求し、その MTO 行を書き出す。

- `__amlo.data` … 直接アクセス、レコード = `(iq, isp)`、中身 `zMLO(1:ldim, 1:ndimMTO)`（complex(8)）
- `__amlo.info` … `ldim, ndimMTO, nqibz, nspx, recl` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)`

`lmf --writeham --mlo`（= `job_mlo` の第 1 段）で生成される。

### 段 1' — `cmlo` の書き出し  【未着手・ここが次の一手】

式 (13) が要るのは `zMLO` ではなく **`cmlo`$_{i\alpha}(q)$**（バンド × MLO、両添字とも $q$ 非依存）。
`Hreduction` の中で作られているが外へ出していない。段 1 と同じ機構で出す。

- `__cmlo.data` … レコード = `(iq, isp)`、中身 `cmlo(1:nband, 1:ndimMTO)`
- `__cmlo.info` … `nband, ndimMTO, nqibz, nspx, recl` / `nskip` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)`

**注意 — バンド索引の整合**: $\Sigma^{\psi}_{ij}(q)$ は `emax_sigm` 以下のバンド、
`cmlo` は MLO の窓（`nskip` 以上）で張られる。**両者の $i$ の原点と範囲を `info` に明記して突き合わせる**。
窓の外のバンドは式 (13) で重み $\bar\theta$ により**連続に**落ちる（鋭い打ち切りではない）。

### 段 2 — $\Sigma^{\rm MLO}(R)$ の生成

入力 `sigm`（または $\Sigma^{\psi}$ を直接）+ `__cmlo.data/.info`、出力 `SigRsMLO`。

1. 式 (13) で $\Sigma^{\rm MLO}(q) = (c^{\rm MLO})^\dagger \Sigma^{\psi}(q)\, c^{\rm MLO}$。**MTO を経由しない。**
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
式 (16) の $A(k) = S^{\rm PMT}(k)\,z^{\rm MLO}(k)$ と $O^{\rm MLO}(k)$ の 2 つだけ。

- $z^{\rm MLO}(k)$ … **内挿しない**。その $k$ で式 (10)–(12) の射影を実行する
  （`Hreduction` は任意 $k$ リストで走る）。$k$ ごとに対角化 1 回分の追加コスト。
- **循環の断ち方**: $\tilde\chi_\alpha(k)$ の定義に使う band 多様体は
  「**直前の反復の $H$ の固有関数**」で固定する。$\Sigma^{\rm MLO}(q)$ を作ったときと同じ処方にすること（§3.2.4 の整合条件）。
  種が固定なのでゲージは自動的に揃う。
- $O^{\rm MLO}(k) = (z^{\rm MLO})^\dagger S^{\rm PMT}(k) z^{\rm MLO}$ を同じ場所で作る。
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
   ずれると式 (13) が静かに壊れる。検証 2（メッシュ点での厳密性）で必ず捕まるので、そこで担保する。
2. **既約 $q$ → 全 BZ の回転** — 段 2-2。`ix(α)` と回転則の対応を誤ると全部壊れる。同じく検証 2 で捕まる。
3. **(B) の循環** — 段 4。$\tilde\chi(k)$ の定義に使う多様体を「直前の反復の $H$」に固定して断つ。
   これが不安定なら (A) のバンド経路だけを採り、SCF は従来の MTO 内挿のまま、という使い分けもあり得る。
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
