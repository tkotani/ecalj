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
hsfp0/hgw : <psi^PMT_i| Sigma |psi^PMT_j>        固有状態基底、既約 q
    |
hqpe_sc   : QSGW の静的 Sigma を組み、固有ベクトルで MTO 基底へ変換
    |       APW 成分は捨てる (mtosigmaonly)
    v
sigm.<sname> : <chi^MTO_mu| Sigma-Vxc |chi^MTO_nu>   ndimsig = nlmto (= ldim)
    |
m_rdsigm2_init : 全 BZ へ展開 (rdsigm2/hamfb3k) -> FFT -> hrr(R)    <- ここが 230^2 x nk^3
    |
bloch2(k)  : WS 最短ベクトルで Bloch 和 -> Sigma(k)                 <- ここで中点がオーバーシュート
    |
getsenex   : PMT ハミルトニアンへ加算
```

## 3. 設計

### 3.0 【2026-09-24 訂正】現状の展開はすでに非直交の双対基底 — 変えるのは部分空間だけ

当初「`mtosigmaonly` が APW を捨てている」と書いたが**誤り**だった。
`getsenex`（`rdsigm2.f90:18`）は

```fortran
ovlmtoi = ovlm(1:ndimsig,1:ndimsig)            ! O_MTO
call matcinv(ndimsig,ovlmtoi)                  ! O_MTO^-1
ovliovl = matmul(ovlmtoi, ovlm(1:ndimsig,1:ndimh))
senex   = matmul(conjg(transpose(ovliovl)), matmul(sene, ovliovl))
```

すなわち

$$\hat\Sigma \;=\; |{\rm PMT}\rangle\,\langle{\rm PMT}|{\rm MTO}\rangle\,O_{\rm MTO}^{-1}\;
\Sigma^{\rm MTO}\;O_{\rm MTO}^{-1}\,\langle{\rm MTO}|{\rm PMT}\rangle\,\langle{\rm PMT}|$$

を既に実装しており、**APW ブロックも埋まる**。`mtosigmaonly` は「Σ の行列要素を MTO 部分空間で保持する」
という意味であって、「APW に効かない」ではない。

したがって **変えられるのは「どの部分空間で Σ を保持し、内挿するか」だけ**である。
MTO（230、非直交、EH/EH2 が準線形従属、裾が長い）→ MLO（154、局在、固定種）に置き換える。
**展開の枠組みは一切変えない**ので、`getsenex` の 4 行を差し替えるだけで済む。

### 3.1 記号

- $\chi^{\rm MTO}_\mu$ … MTO 基底（$\mu = 1\ldots L$、$L$ = `ldim` = 230）
- $|w_i\rangle = \sum_\mu \chi^{\rm MTO}_\mu A_{\mu i}$ … MLO（$i = 1\ldots M$、$M$ = `ndimMTO` = 154）。
  $A(q)$ = `zMLO(1:ldim, 1:ndimMTO)`、`Hreduction` が返す。**段 1 で `__amlo.data` に書き出し済み**。
- $\Sigma^{\rm MTO}(q)$ … `sigm` の中身

### 3.2 MLO 表現

$$\Sigma^{\rm MLO}(q) = A^\dagger(q)\,\Sigma^{\rm MTO}(q)\,A(q) \qquad (M \times M)$$

計量は不要（`sigm` は $\langle\chi|\Sigma-V_{xc}|\chi\rangle$ そのもの）。自由度は $L^2 \to M^2$ で **45 %**。

### 3.3 内挿

$\Sigma^{\rm MLO}(q)$ を FFT して $\Sigma^{\rm MLO}(R)$ とし、任意 $k$ では
`HamRsMLO` と同じ WS 最短ベクトル・重みで Bloch 和を取る。
**内挿されるのはここだけ**。

### 3.4 `getsenex` の差し替え

$$\hat\Sigma = |{\rm PMT}\rangle\,\langle{\rm PMT}|w\rangle\,O_{\rm MLO}^{-1}\;
\Sigma^{\rm MLO}(k)\;O_{\rm MLO}^{-1}\,\langle w|{\rm PMT}\rangle\,\langle{\rm PMT}|$$

コード上は `ovlm(1:ndimsig,1:ndimh)` を $A^\dagger(k)\,\mathrm{ovlm}(1{:}L,1{:}n_{\rm dimh})$ に、
`ovlm(1:ndimsig,1:ndimsig)` を $A^\dagger(k)\,\mathrm{ovlm}(1{:}L,1{:}L)\,A(k)$ に置き換えるだけ。
$M \times M$ の逆行列（154³）は無視できるコスト。

### 3.5 $A(k)$ をどう作るか

$A(q)$ は Σ メッシュ点でしか無いので、任意 $k$ には

$$C(R) = \mathrm{FFT}[A(q)], \qquad A(k) = \sum_R C(R)\,e^{ikR}$$

**ゲージは projection gauge** — MLO は固定した MTO 種 $F^{\rm MTO}_k$ を band 多様体へ射影して作る
（`m_hreduction.f90:317`）ので、固有ベクトルの任意位相は入らない。したがって $C(R)$ は意味を持つ。
**要検証**（§5 の検証 1・2）。

反復の中では **直前の反復の $A$ を使う**（user 2026-09-24）。QSGW が $\Sigma$ 自体についてやっていることと
同じで、自己無撞着に至れば整合する。

### 3.6 何を失うか

MLO 部分空間の外の $\Sigma$。LiTi₂O₄ の 154 軌道模型では

| 帯 | メッシュ点での再現 |
|---|---|
| O 2s / O 2p / t2g / eg（E_F+4 eV まで） | **0.1〜0.5 meV** |
| 格子間バンド（E_F+4 eV 以上） | 33〜306 meV |

### 3.7 PAW チャネル案（将来）

MLO の代わりに MT 球内の部分波 φ, φ̇ のチャネル（`@MNLA_CPHI`、`ndima` = 700）へ射影する案もある。
利点は $B(k)$（augmentation 係数）が **lmf が任意 $k$ で厳密に作れる**ため $A(k)$ の内挿が要らないこと。
欠点は次元が 700 と多いこと（l ≤ 3 で 448、l ≤ 2 で 252、原子別に切れば 〜150）。
**MLO 案で §5 の検証 1・2 が通らなかった場合の代替**として残す。

## 4. 実装手順

### 段 1 — $A(q)$ の書き出し  【2026-09-24 完了 — 実装済みはここまで】

`m_HamPMT.f90` の `HreductionIqibz` ループで `zMLO` を常時要求し、その MTO 行を書き出す。

- `__amlo.data` … 直接アクセス、レコード = `(iq, isp)`、中身 `zMLO(1:ldim, 1:ndimMTO)`（complex(8)）
- `__amlo.info` … `ldim, ndimMTO, nqibz, nspx, recl` / `ix(1:ndimMTO)` / `qibz(1:3,1:nqibz)`

既存の SOC 用 `zMLO` をそのまま使うので追加コストはほぼゼロ。
`lmf --writeham --mlo`（= `job_mlo` の第 1 段）で生成される。

### 段 2 — $\Sigma^{\rm PAW}(R)$ の生成

新規ツール（または `hqpe_sc` の後段）。入力 `sigm`, `__amlo.data/.info`、出力 `SigRsMLO`。

1. `hqpe_sc` で $\Sigma^{\rm PAW}_{ab}(q) = \sum_{ij}\overline{c_{ai}}\,\Sigma_{ij}\,c_{bj}$ を作り、
   `sigm` の代わりに（または併記して）書き出す。`cphi` は GW 側に既にある。
2. 既約 q → 全 BZ へ展開。**$P_a$ は原子中心の球面調和なので回転則が明快**
   （MTO/MLO のように「行と列で回転則が違う」問題が無い）。
3. FFT → $\Sigma^{\rm PAW}(R)$、WS 最短ベクトルの対リストで保存。

### 段 3 — `bloch2` の差し替え

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
