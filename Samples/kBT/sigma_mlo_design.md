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

## 2. 現状のデータフロー

```
hsfp0/hgw : Sigma^psi_ij(q) = <psi^PMT_i| Sigma |psi^PMT_j>     固有関数基底、既約 q      … 式 (7)
    |
hqpe_sc   : QSGW の静的 Sigma を組み、MTO 基底の行列へ変換                              … 式 (9)
    v
sigm.<sname> : Sigma^MTO_munu(q)     ndimsig = nlmto (= L = 230)
    |
rdsigm2   : 対称操作で全 BZ へ展開 (hamfb3k)
    |
fftz3     : 実空間へ    Sigma^MTO_munu(R)   =  hrr                 <- 自由度 L^2 x nk^3
    |
bloch2(k) : WS 最短ベクトル・重み平均で任意 k へ Bloch 和                               … 式 (1)
    |
getsenex  : 双対展開で PMT 基底へ広げ、H に加える                                        … 式 (2)-(4)
```

**内挿されるのは式 (9) の行列要素**であり、$\mu\nu$ 成分ごとに独立にフーリエ内挿される:

$$\Sigma^{\rm MTO}_{\mu\nu}(R) = \frac{1}{N_q}\sum_{q} \Sigma^{\rm MTO}_{\mu\nu}(q)\,e^{-iqR},
\qquad
\Sigma^{\rm MTO}_{\mu\nu}(k) = \sum_{R} \Sigma^{\rm MTO}_{\mu\nu}(R)\,\overline{e^{ikR}}\tag{1}
$$

（上線は縮退した最短ベクトル群についての平均、`m_shortn3.f90:gennlat` の `nqwgt`。）

**問題はここ**。6³ なら既約 q は 16 点しかないのに、式 (1) は $L^2 = 230^2$ 個の要素それぞれについて
$n_k^3 = 216$ セル分の実空間自由度を持たせている。しかも MTO 基底は非直交で EH/EH2 が準線形従属、
包絡関数の裾も長い。結果として、**メッシュ点の間で線形独立性の低い自由度がふらつく**（§1 の表）。

## 3. 設計

### 3.0 記号 — 基底・双対基底・MLO・固有関数を区別する

| 記号 | 意味 | 次元 / 索引 |
|---|---|---|
| $\chi^{\rm PMT}_m$ | **PMT 基底関数**（MTO + APW） | $m = 1\ldots n_{\rm dimh}$、**APW の本数は $k$ に依存** |
| $\chi^{\rm MTO}_\mu$ | その MTO 部分（**下付き = 共変**） | $\mu = 1\ldots L$（`ldim` = 230）、**$k$ 非依存** |
| $\chi_{\rm MTO}^{\mu}$ | その**双対（反変）基底**、$\langle\chi_{\rm MTO}^{\mu}|\chi^{\rm MTO}_\nu\rangle=\delta^\mu_\nu$ | 上付きで区別 |
| $\chi^{\rm MTO}_{ix(k)}$ | **MLO の種**。`ix(k)` 番目の MTO 基底関数（固定） | 新しい記号は要らない |
| $\tilde\chi_k = P\,\chi^{\rm MTO}_{ix(k)}$ | **MLO**（コードの `F^MLO`）。band 多様体への射影 | $k = 1\ldots M$（`ndimMTO` = 154） |
| $\psi^{\rm PMT}_i$ | **固有関数**（基底ではない） | |
| `cmlo`$_{ik} = \langle \psi^{\rm PMT}_i|\tilde\chi_k\rangle$ | MLO の**固有関数**基底での係数 | |
| `zMLO`$_{mk}$ | MLO の **PMT 基底**での係数、$\tilde\chi_k=\sum_m \chi^{\rm PMT}_m z^{\rm MLO}_{mk}$ | $(n_{\rm dimh}, M)$ |
| $P_a$ | **PAW（augmentation）チャネル** φ, φ̇ | $a = 1\ldots N_a$（`ndima` = 700）、**$k$ 非依存** |
| $B_{am}(k)$ | PMT 基底関数の augmentation 係数 | lmf が任意 $k$ で厳密に作る |

コードの `F^MTO_k` は $\chi^{\rm MTO}_{ix(k)}$、`F^MLO_k` は $\tilde\chi_k$ に対応する。

$\Sigma$ は `sigm` に $\Sigma^{\rm MTO}_{\mu\nu}(q)=\langle\chi^{\rm MTO}_\mu|\hat\Sigma|\chi^{\rm MTO}_\nu\rangle$ として入っている（$L\times L$）。

### 3.1 現状の展開はすでに非直交の双対基底【2026-09-24 訂正】

当初「`mtosigmaonly` が APW を捨てている」と書いたが**誤り**。`getsenex`（`rdsigm2.f90:18`）は

```fortran
ovlmtoi = ovlm(1:ndimsig,1:ndimsig)            ! O_MTO
call matcinv(ndimsig,ovlmtoi)                  ! O_MTO^-1
ovliovl = matmul(ovlmtoi, ovlm(1:ndimsig,1:ndimh))
senex   = matmul(conjg(transpose(ovliovl)), matmul(sene, ovliovl))
```

を計算している。**すべて基底関数**（波動関数ではない）の話であることに注意。記号を決める:

- $O^{\rm MTO}_{\mu\nu}(k) = \langle\chi^{\rm MTO}_\mu|\chi^{\rm MTO}_\nu\rangle$ … MTO 同士の重なり（`ovlm(1:L,1:L)`）
- $S_{m\mu}(k) = \langle\chi^{\rm PMT}_m|\chi^{\rm MTO}_\mu\rangle$ … PMT 基底と MTO 基底の重なり（`ovlm(1:L,1:ndimh)` の転置）

MTO 基底は非直交なので、**双対（反変）基底**

$$|\chi_{\rm MTO}^{\mu}\rangle \;=\; \sum_\nu |\chi^{\rm MTO}_\nu\rangle\,(O^{\rm MTO})^{-1}_{\nu\mu},
\qquad \langle\chi_{\rm MTO}^{\mu}|\chi^{\rm MTO}_\nu\rangle = \delta^{\mu}_{\nu}\tag{2}
$$

を使って演算子を組む:

$$\hat\Sigma \;=\; \sum_{\mu\nu} |\chi_{\rm MTO}^{\mu}\rangle\;\Sigma^{\rm MTO}_{\mu\nu}\;\langle\chi_{\rm MTO}^{\nu}|\tag{3}
$$

その **PMT 基底での行列要素**が `senex`:

$$\big[\hat\Sigma\big]_{mn} \;=\; \langle\chi^{\rm PMT}_m|\hat\Sigma|\chi^{\rm PMT}_n\rangle
\;=\; \sum_{\mu\nu} T^{*}_{\mu m}\;\Sigma^{\rm MTO}_{\mu\nu}\;T_{\nu n},
\qquad T_{\mu n} \;=\; \sum_\rho (O^{\rm MTO})^{-1}_{\mu\rho}\,S^{*}_{n\rho}\tag{4}
$$

（$T$ = コードの `ovliovl`。）$m,n$ は **MTO と APW の両方**を走るので、**APW ブロックも埋まる**。
`mtosigmaonly` は「$\Sigma$ の行列要素を MTO 部分空間で保持する」の意であって、
「APW に効かない」ではない。

→ **変えられるのは「どの部分空間で $\Sigma$ を保持し、内挿するか」だけ**。展開の枠組みは変えない。

### 3.2 MLO を部分空間に使う場合【障害あり】

$\tilde\chi$ で保持するなら

$$\Sigma^{\rm MLO}_{kl}(q) = \langle \tilde\chi_k|\hat\Sigma|\tilde\chi_l\rangle
= \big(D^\dagger \Sigma^{\rm MTO} D\big)_{kl},\qquad
D(q) = O_{\rm MTO}^{-1}\,\langle{\rm MTO}|{\rm PMT}\rangle\, z^{\rm MLO}\tag{5}
$$

戻すときは

$$\hat\Sigma = |{\rm PMT}\rangle\,S_{\rm PMT} z^{\rm MLO}\,O_{\rm MLO}^{-1}\;\Sigma^{\rm MLO}(k)\;
O_{\rm MLO}^{-1}\,(z^{\rm MLO})^\dagger S_{\rm PMT}\,\langle{\rm PMT}|,
\qquad O_{\rm MLO} = (z^{\rm MLO})^\dagger S_{\rm PMT} z^{\rm MLO}\tag{6}
$$

**障害**: 任意 $k$ で $z^{\rm MLO}(k)$ が要るが、$z^{\rm MLO}$ の行は PMT 基底で、
**APW の本数が $k$ ごとに違う**（$|k+G| <$ cutoff の $G$ 集合が変わる）。
したがって **$C(R) = \mathrm{FFT}[z^{\rm MLO}(q)]$ が定義できない**。
MTO 行だけなら $k$ 非依存で FFT できるが、それは MLO の一部でしかない
（切り捨てると $\tilde\chi$ ではない別物になる）。

→ 式 (5)(6) の道は、この障害を回避しない限り成立しない。

回避案（いずれも未検証）:
1. MLO の APW 成分を無視し、MTO 成分だけで近似する。$\tilde\chi$ の APW 重みが小さい系でのみ可。
2. **MLO を PAW チャネルで展開する**（→ §3.3）。PAW チャネルは $k$ 非依存なので係数が固定できる。

### 3.3 PAW（augmentation）チャネルを部分空間に使う【本命】

**出発点**。GW（`hsfp0`/`hgw`）が各既約 $q$ で出すのは、**固有関数で挟んだ行列要素**

$$\Sigma^{\psi}_{ij}(q) \;=\; \langle \psi^{\rm PMT}_{iq} \,|\, \hat\Sigma \,|\, \psi^{\rm PMT}_{jq}\rangle\tag{7}
$$

である（式 (7)。QSGW ではこれをエルミート化した静的 $\Sigma$）。演算子としては

$$\hat\Sigma(q) \;=\; \sum_{ij} |\psi^{\rm PMT}_{iq}\rangle\;\Sigma^{\psi}_{ij}(q)\;\langle \psi^{\rm PMT}_{jq}|\tag{8}
$$

固有関数は正規直交（$\langle\psi_i|\psi_j\rangle = \delta_{ij}$）なので、ここには逆行列が要らない。

現状の `hqpe_sc` は、これを **MTO 基底**の行列へ変換している:

$$\Sigma^{\rm MTO}_{\mu\nu}(q) \;=\; \langle\chi^{\rm MTO}_\mu|\hat\Sigma(q)|\chi^{\rm MTO}_\nu\rangle
\;=\; \sum_{ij} \langle\chi^{\rm MTO}_\mu|\psi_i\rangle\,\Sigma^{\psi}_{ij}\,\langle\psi_j|\chi^{\rm MTO}_\nu\rangle\tag{9}
$$

式 (9) が `sigm` の中身で、**途中産物**にすぎない。
$\hat\Sigma$ を別の部分空間で保持したいなら、**MTO を経由せず $\Sigma^{\psi}$ から直接**取ればよい。

$P_a$（MT 球内の φ, φ̇）は **$k$ 非依存の固定索引**であり、$B(k)$ は **lmf が任意 $k$ で厳密に作る**。
したがって §3.2（式 (5)(6)）の障害（APW の本数が $k$ 依存）が最初から無い。

**(1) GW の出力から直接**

固有関数の augmentation 係数 `cphi` は $\psi_i = \sum_a P_a\,c_{ai}$（球内）なので
$\langle P_a|\psi_i\rangle = (\Pi c)_{ai}$（$\Pi$ = 球内の重なり `ppj`）。演算子 $\hat\Sigma = \sum_{ij}|\psi_i\rangle\Sigma^\psi_{ij}\langle\psi_j|$ の PAW 行列は

$$\boxed{\ \Sigma^{\rm PAW}(q) \;=\; (\Pi\,c)\;\Sigma^{\psi}(q)\;(\Pi\,c)^\dagger\ }\tag{10}
$$

`hqpe_sc` が今 式 (9) を作っているところを、式 (10) に差し替える。
`cphi` も `ppj` も GW 側に既にある。

**(2) 内挿**

$\Sigma^{\rm PAW}(q)$ を FFT して $\Sigma^{\rm PAW}(R)$、任意 $k$ で Bloch 和。**内挿されるのはここだけ。**

**(3) 戻す**

$\langle\chi^{\rm PMT}_m|P_a\rangle = (B^\dagger \Pi)_{ma}$ なので、双対展開の $\Pi^{-1}$ が両側で約分して

$$\boxed{\ \Sigma^{\rm PMT}_{mn}(k) \;=\; \big(B^\dagger(k)\,\Sigma^{\rm PAW}(k)\,B(k)\big)_{mn}\ }\tag{11}
$$

式 (11) の **$B(k)$ は厳密**、逆行列も不要。`getsenex` は 1 行になる。

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
| 1 | **$\Sigma^{\rm PAW}(R)$ の減衰** … $\|\Sigma(R)\|$ vs $\|R\|$ をチャネル別に | MTO 表現（09-24 の測定）より速く減衰。BvK セル端で頭打ちしない |
| 2 | **次元の切り詰め** … l の上限を 4/3/2 と変えて $\Sigma$ の値と荒れを見る | 値が動かない範囲で最小の l |
| 3 | **メッシュ点での厳密性** … §3.5 の $\Sigma^{\rm PMT}(k)$ を q メッシュ点で元と比較 | バンドで 1 meV 以内 |
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
