# kBT — 有限温度 QSGW (`tetrakbt` + `t_sigmakbt`) のサンプル

LiTi₂O₄ (金属スピネル) の QSGW を、**電子温度 2000 K** を χ₀ 側と Σ 側の
両方に入れて回したもの。手法の説明は
[ecaljdoc: 有限温度の自己エネルギー計算](https://ecalj.github.io/manual/kBT) を見よ。

> **これは入力と結果を置いてあるだけのサンプルである。**
> GW を 11 反復するので計算が重く、`testecalj` のターゲットにはしていない。
> 手法が何を変えるかを、収束した結果そのもので確認するためのもの。

ブランチ: `t_sigmakbt` (Σ 側の有限温度。χ₀ 側の `tetrakbt` は `main` にある)

---

## 1. 何を計算したか

| 項目 | 値 |
|---|---|
| 系 | LiTi₂O₄ (spinel, 金属) |
| GW の k メッシュ `n1n2n3` | 6×6×6 |
| `deltaq_scale` | 0.1 |
| `tetrakbt` / `t_tetrakbt` | `true` / 2000.0 K (χ₀ = $W$ 側) |
| `t_sigmakbt` | 2000.0 K ($\Sigma$ 側、`t_tetrakbt` と同じ温度) |
| QSGW 反復数 | 11 |
| $k_BT$ | 0.01267 Ry = 0.172 eV |

`ctrlg.liti2o4.toml` の該当部分:

```toml
[gw]
n1n2n3       = [6, 6, 6]
deltaq_scale = 0.1
tetrakbt     = true      # chi0 を有限温度テトラヘドロンで (method B')
t_tetrakbt   = 2000.0    # 2000 K
t_sigmakbt   = 2000.0    # Sigma 側の FD 有限温度 = t_tetrakbt
esmr         = 0.01      # 数値的な極の平滑化 (有限温度は t_sigmakbt が担う)
```

**`t_sigmakbt` は `tetrakbt = true` と必ず併用すること。**
Σ 側が読む `EFERMI_kbt` を書くのは `tetrakbt` を有効にした `heftet` なので、
単独では警告が出て有限温度が適用されない。

---

## 2. ディレクトリ

```
input/     ctrlg.liti2o4.toml   入力 (上記の [gw] を含む)
           PB.liti2o4.toml      product basis
           syml.liti2o4         バンド線 (Γ-X-U|K-Γ-L-W-X)
           rst.liti2o4.lda      LDA 収束済みの rst (ここから gwsc を始めた)

results/   sigm.liti2o4         11 反復後の収束した QSGW 自己エネルギー (13 MB)
           EFERMI, EFERMI_kbt   T=0 と有限温度の Fermi 準位
           finiteT_evidence.txt 実行ログからの抜粋 (下記 §4)
           QPU, QPU.1run .. QPU.11run    各反復の QP エネルギー
           band_iterations.png  11 反復を 1 枚に重ねた図
           band_iter/band_iter1.png .. 11.png   反復ごとのバンド
           bnd_iterations.tar.gz  上の図の元データ (展開すると
                                  bnd/iter1 .. iter11/bnd00[1-6].spin1 +
                                  bandplot.isp1.glt。展開後 27 MB)

plots/     T_dependence.png          T = 0/1000/2000/3000 K の比較 (別 run 群、dq=0.3)
           deltaq_artifact_3000K.png 3000 K の異常が deltaq 由来である証拠

README_kt1_original.md   計算を回した作業ディレクトリ (kt1) の元 README
```

---

## 3. 結果

### 3.1 QSGW 反復の収束

![band iterations](results/band_iterations.png)

紫 (反復 1) から黄 (反復 11) へ。反復 3–4 以降は $E_F$ 近傍でも線が重なり、
金属にもかかわらず QSGW が振動せずに収束している。これが `tetrakbt` +
`t_sigmakbt` を入れる目的である。反復ごとの個別の図は
`results/band_iter/band_iter<N>.png`。描き直すなら:

```bash
cd results && tar xzf bnd_iterations.tar.gz   # -> bnd/iter1 .. iter11/
```

### 3.2 温度依存 — 何が直るか

![T 依存](plots/T_dependence.png)

T = 0 / 1000 / 2000 / 3000 K を重ねたもの (QSGW 第 1 反復、6³、`deltaq_scale = 0.3`。
共通の静電ゼロで揃え、T=0 の E_F を原点に取った)。これは本サンプルとは
**別の run 群** (`deltaq_scale` が 0.3) からのもの。

T=0 (`tetrakbt` off = 従来の `lindtet6`) では Ti-3d t2g 帯が Γ–L 上で
**−3.9 eV まで落ち込む**。これは分散ではなく q→0 head の破綻で、温度を
上げると単調に消える:

| T | Γ–L 上 (x≈2.7) の t2g 帯の底 |
|---|---|
| 0 (`lindtet6`) | −3.91 eV |
| 1000 K | −1.31 eV |
| 2000 K | −0.66 eV |
| 3000 K | −0.50 eV |

これが有限温度を入れる理由である。なお method B′ が T→0 で `lindtet6` に
厳密に戻ること自体は変わらない (T=0 で壊れているのはテトラヘドロン法ではなく、
その上での金属の q→0 の扱いのほう)。

> 500 K の run もアーカイブにあるが、**そのバンドプロットは E_F を誤って拾っており
> (band ステップの `efermi.lmf` が 0.158 Ry、他の run は 0.27–0.28 Ry)、
> 分散自体も揃わない**ので上の図からは外してある。

### 3.3 3000 K の K–Γ 異常は `deltaq_scale` のアーティファクト

![deltaq artifact](plots/deltaq_artifact_3000K.png)

上の図の 3000 K は Γ–L では収まっているのに K–Γ の中央に別のスパイクが残る。
`deltaq_scale = 0.3` (赤) では **−1.24 eV のスパイク**が出るが、
`deltaq_scale = 0.1` (青) にすると**消える**。

これは `tetrakbt` / method B′ 本体の問題ではなく、**offset-Gamma の q→0
head が、高温 × 大きい `deltaq` で破綻する**もの。高温を使うときは
`deltaq_scale` を小さく取ること。本サンプルが `deltaq_scale = 0.1` なのはこのため。

---

## 4. 有効になっていることの確認

`results/finiteT_evidence.txt`:

```
--- chi0 側 (tetrakbt) ---
 tetrakbt_init: T[K], kbt[Ry], kbt[eV] 2000  0.12667E-01  0.17234E+00

--- Sigma 側 (t_sigmakbt) ---
 sigmakbt_setup: t_sigmakbt[K]=   2000.0 kbt[Ry]=  0.12667E-01 ef<-EFERMI_kbt=    0.267452

--- EFERMI (T=0) / EFERMI_kbt (有限温度) ---
  0.277063046836211D+00 ... ! efermi bandgap are obtained by heftet
  0.267452126929387D+00 ... ! efermi_kbt (finite-T, tetra-NOS conv) by heftet
```

2 行目の `ef<-EFERMI_kbt` が、Σ 側が **T=0 の `EFERMI` (0.2771 Ry) ではなく
有限温度の `EFERMI_kbt` (0.2675 Ry) を使っている**ことを示す。
この 2 つが食い違ったままだと χ₀ と Σ が別の Fermi 準位を見ることになる。

---

## 5. 自分で回すなら

```bash
mkdir work && cd work
cp ../input/ctrlg.liti2o4.toml ../input/PB.liti2o4.toml ../input/syml.liti2o4 .
cp ../input/rst.liti2o4.lda rst.liti2o4
lmfa liti2o4
gwsc 11 -np <NP> liti2o4      # 重い。kt1 で 60 プロセス
```

バンドだけ見たいなら GW をやり直さず、収束した自己エネルギーを使えばよい:

```bash
cp ../results/sigm.liti2o4 sigm.liti2o4
job_band liti2o4 -np <NP> -vnspin=1
```
