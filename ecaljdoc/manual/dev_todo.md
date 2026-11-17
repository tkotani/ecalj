# TODO（開発の残り） {#dev-todo}

::: info このページについて
ecalj の開発で、まだ済んでいないことの一覧です（2026-11-17 に作った）。
正本はリポジトリの [`MD/TODOandQuestion.md`](https://github.com/tkotani/ecalj/blob/main/MD/TODOandQuestion.md) で、
このページはそこから、外の人にも役に立つものを選んで写しました。計算機ごとの作業の進み具合は写していません。食い違うときは正本が正しい。
各項目の括弧の日付は、見つけた日です。印: **【未着手】**、**【やりかけ】** 手を付けて途中、**【判断待ち】** 方針を決めていない、**【残す】** 急がない。
:::

## 各ページにある未解決の項目 {#in-pages}

ここには写さず、各ページの節を指す。

- MLO: [残る問題](./mlo#mlo-problems)（固定の $w$ で半導体と金属を同時に満たせない、など）
- MLO-gwsc: [MLO-gwsc](./mlo_gwsc) の節「未解決・今後」（`pwmode = 11` との比べ、反強磁性の対称性との併用、MLO の本数の影響、など）
- kBT: [kBT](./kBT) の節「残っている課題」（金属での `t_sigmaw` の選び方、絶縁体の `EFERMI_kbt`、Bose 項、など 8 項目）
- 書きかけの理論のページ: [光学的性質](../theory/optical_properties)、[Basic II](../theory/basic2)、[GW で用いる波動関数の表現](../theory/gw_basis)、
  [積基底](../theory/gw_product_basis)、[GW で用いる行列要素](../theory/gw_zmel)

## MLO {#mlo}

- **【未着手】EuO に EH2 を足すと模型を作る所で止まる**（2026-10-01）。EH と EH2 の両方をシードにすると、`zhev_tk4` が一次従属に近い向きを落とす。
  4f の原子は [模型の選び方](./mlo#mlo-choice) の基準 2 から除く決まりなので急がない
- **【残す】マグノンの Goldstone の条件の倍率 η ≠ 1 の原因**（2026-10-02）。Löwdin の MLO で Fe 1.24、FeCo 1.28、Ni 1.75。Ni が大きい理由は未確認
- **【未着手】Fe のマグノン: 小さい q で 2 割高い**（2026-10-02）。Löwdin の MLO は q ≤ 0.3 で Wannier 版より 2 割高い（q = 0.1 で 0.085 対 0.068 eV）。
  q ≥ 0.4 は合う。実験（中性子散乱）の分散との比較もまだ（データはリポジトリに無い）。
  サンプルは [`Samples/Magnon/Fe_mlo_magnon`](https://github.com/tkotani/ecalj/tree/main/Samples/Magnon/Fe_mlo_magnon)
- **【残す】空格子球（ES）の自動化**（2026-10-01）。[`ctrlg_addes.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/ctrlg_addes.py) を ctrlg の生成に組み込むか。
  板（厚い真空のある構造）では、表面の原子の外に積層を延ばして ES（s,p）を置くと MLO が通る（2026-10-08 に規則にした）。この置き方はまだ `ctrlg_addes.py` に無い
- **【未着手】ES で対称性が下がる**（2026-10-05）。`ctrlg_addes.py` は重ならない空隙を大きい順に採るので、対称性で等価な空隙の一部だけを採ることがある。
  既約な q 点が増え、GW が 2〜6 倍重くなる。等価な空隙はまとめて採るようにする
- **【未着手】`gwinit` が ES（z = 0）の `mlo_lm` を基底に合わせて書く**（2026-10-05）。いまは s,p を書くので、s だけの ES では基底に無い p が入る

## 対称性と GW {#gw}

- **【未着手】EIBZ を復活させ、Σ（と χ0）を q を保つ操作の群で対称化する**（2026-10-02）。
  対称操作の並べ方を変えると、AF NiO（R-3m、操作 12）の QSGW の 1 反復目の QP が最大 15 meV 変わる。帯の数の上限や `emax_sigm` の切れ目のような、
  回転と両立しない切り捨てがあるためで、バグではない（試験の参照は今の並べ方の値にした）。Σ(q) を $\Sigma\leftarrow\frac{1}{N_g}\sum_g U_g\Sigma U_g^\dagger$ で平均すれば
  この差が消え、k の和も減る。EIBZ の [`eibzgen.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/eibzgen.f90) は残っている。
  まず numpy で NiO 4³ の Σ を対称化して、差が消えるかと時間を見る
- **【未着手】時間反転の無い GW（`npm = 2`）**（2026-10-02）。スピン軌道を GW の中に入れた磁性体で空間反転も無いと、負の振動数の χ0・W・Σc が別に要る。
  いまは χ0 が `x0kf_zxq` の冒頭で、Σ が `m_sxcf_sc_count.f90` の `npm=2 need to be examined` で止まる。
  手順の案: (a) 空間反転があれば、負の振動数の $W(\mathbf q)$ を $W_{\mathbf G\mathbf G'}(\mathbf q,-\omega)=W_{-\mathbf G',-\mathbf G}(-\mathbf q,\omega)$ で作る、
  (b) 無ければ振動数の積分を負の側まで通す。オルタ磁性体（RuO₂、MnTe は空間反転があるので (a)）で要る
- **【未着手】コアの交換項を積基底（MPB）を使わずに計算する**（2026-11-15）。いまは `gwsc` の反復ごとにコア用の MPB を作ってから Σ^Cx を出している
  （[実装メモ](./dev_impl) の (I15) の下）。コアは MT の中だけにあるので、MT の中の積分として直接計算できるはず
- **【未着手】`hgw` の残り**（2026-04）。ノードの中で W を共有する（`MPI_Win_allocate_shared`）、`hsfp0_sc` の Sx・コアの交換も `hgw` に入れる（優先度は低い）
- **【未着手】`lmf --jobgw=1` のメモリ**（2026-10-01）。`GEIGpart` の `ppovl(ngp,ngp)` を各ランクで持つ（32·ngp² バイト、胞の体積の 2 乗）。
  案: ngp から並列数を決める、Cholesky の因子だけ持つ、$O\cdot x$ を FFT で作り反復法で解く

## 有限温度 {#finite-t}

- **【未着手】有限温度の四面体法と T = 0 の四面体法の関係を確かめる**（2026-10-02）。`t_tetrakbt > 0`（`m_tetrakbt`）の重みが、T → 0 で
  `t_tetrakbt = 0`（`tetwt5`）の重みに一致するか。同じ k メッシュで χ0 の虚部と QP のエネルギーを、T = 300、100、30、10 K と下げて比べる。
  Σ 側の幅 `t_sigmaw` との組み合わせも。ほかの課題は [kBT](./kBT) の「残っている課題」

## 入力、k 点、SCF {#input}

- **【やりかけ】k 点のメッシュを間隔で決める規則**（2026-10-05）。いまは `ctrlgenToml.py` の既定が 8×8×8、GW が 4×4×4 の一律。
  点の間隔 Δk で $n_i=\max(n_{\min},\lceil|\mathbf b_i|/\Delta k\rceil)$ とし、lmf は Δk = 0.26 Å⁻¹（Si で 8）、GW は 0.34 Å⁻¹（Si で 6、8 原子の胞で多くは 4）と決めた。
  `kmesh.py` はブランチ `kmesh` にあり、main にまだ入れていない。`ctrlgenToml.py` が `[bz] nkabc` と `[gw] n1n2n3` の既定をこれで書くようにする
- **【未着手】金属の GW の k 点の規則**（2026-10-05）。Na・Al（sp）と Fe・Ni・Cu（3d）で GW のメッシュを振り、QP の準位・d バンドの幅・交換分裂が
  0.05 eV に収まる段を見る。$k_F$ で整理できるか、d の状態密度で整理できるかで、金属だけ段を上げる規則を決める
- **【未着手】混合のパラメータ b の自動の調整を lmf の中で完結させる**（2026-10-05）。いまは Python（[`pylib/dft.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/dft.py) の `bmix_reduction`）が、
  lmf が落ちるたびに b を下げて lmf を初めからやり直す。遅く、lmf を直接回す人には効かない。
  案: `m_lmfp.f90` の SCF のループで発散の兆候（qdiff が続けて増える、NaN）を見たら、反復の初めの密度に戻し、混合の履歴を捨て、b を半分にして続ける。
  save ファイルの行末に `b=0.20->0.10 at it 7` のように書く
- **【判断待ち】分子の設定**（2026-11-14）。`ctrlgenToml.py --systype=molecule` は `tetra = true`、`w = 0.01` を書くが、
  `toml_comments.py` のコメントと [Samples/AtomDimer/N2](https://github.com/tkotani/ecalj/tree/main/Samples/AtomDimer) は `tetra = false`、`w = 0.001`。
  分子を計算するときに両方を比べて、どちらかに合わせる
- **【判断待ち】試験の入力の温度**。試験は `t_tetrakbt = 262`、`t_sigmaw = 0`。テンプレート（300、300）に揃えると、gas_gwsc・fe_gwsc などの参照が動く
- **【未着手】`Samples/MLOsamples/RuO2` の `rst`・`dmats`** は、`pwmode = 11` の LDA+U に誤りがあった時期（2026-03-30〜09-30）に作ったもの。作り直す

## コード {#code}

- **【やりかけ】PMT の静電・コアまわりのコメント**（2026-11-13）。誤ったコメントは直した。残り: (a) コメントの式番号は論文（JPSJ 84, 034702）のものなので、
  [ecalj のバンド計算法](../theory/pmt) の式番号も併記する、(b) `locpot` の後ろの、コメントにした古い `locpt2`（約 200 行）を外す
- **【未着手】異常ホール伝導度（AHC）のコードの取り込み**（2026-10-02）。前の AHC（`hahc`、`job_AHC`）は Wannier の経路と一緒に外した（git のタグ `last-wannier`）。
  新しいコードを入れ直す
- **【未着手】AFTEST（反強磁性のモーメントの固定）の残り**（2026-10-01）。(a) モーメントを下げる向きで更新が行き過ぎる、(b) 対がサイト 1・2、ブロック 1・2 の決め打ち、
  (c) lmf の起動のたびに場を更新して `mmagfield.aftest` を書き直す（`job_band` でも）
- **【未着手】`SRC/` に残っている使われていないらしいファイル**（2026-10-02）。`SRC/.#Memo4rotation`、`SRC/exec/` の `cmake_install.cmake`、`hello.py`、`platform`、
  `lmf2.py`、`lmchk.py`、`pylmfa`、`pysample`、`ohtaka`、`epsPPd`、`epsPPsaito`、`job_senefbz`、`readeps_dig2.py`、`auto_kauto.py`。使われていないことを確かめて外す
- **【未着手】コードの中の確かめの印**: 力の式の項に `need check`（`m_bndfp.f90` の `frcbandsym`、`fh_rv`）、
  `m_bandcal.f90` の「m_igv2x を、どのランクからもどの k 点のデータも読めるように直す」

## GW1500 と Materials Project {#gw1500}

GW1500 は Materials Project の約 1500 物質の QSGW のデータベース（[DOSnpSupplement](https://github.com/tkotani/DOSnpSupplement)）。

- **【未着手】`auto_mpquery.py` が今の Materials Project で動かない**（2026-10-01）。新しい ID の形で `mp_api` の検証が止まる。`mp_api` を上げるか、REST を直接読む
- **【判断待ち】構造の疑わしい物質**（2026-10-01）。副格子を抜き出した MP の項目が選ばれていた。`INVALID_STRUCTURE`（12）・`SUSPECT_STRUCTURE`（4）を集合から外すか（いまは注記だけ）。
  選ぶときに `auto_mpquery.py` で弾くか印を付ける
- **【未着手】8³ のメッシュがバンドの極値を取り逃がす物質**（51 物質）。黒鉛型の炭素は半金属なのにギャップが出る。ギャップを経路の点も含めて決めるか、メッシュを細かくする
- **【未着手】QSGW80 が LDA より小さい物質の確かめ**: Mo₂O₆（mp-796276）、CdSnAs₂（mp-3829、実験 0.26 eV）
- **【判断待ち】板の物質の QSGW80 は真空の厚さに依る**（2026-10-08）。同じ PbS の板で LDA 1.728 / 1.738 eV、QSGW80 3.725 / 4.152 eV（真空 13.7 / 19.4 Å、周期的な像からの遮蔽）。
  データベースに注記するか、真空を揃えて回し直すか
- **【未着手】MLO が通らない物質の見直し**（2026-10-06）。データベースの「Later list」の物質を一件ずつ見る。C₃N₄ mp-570572、Rb₂Sc₂O₄ mp-7650、H₂ mp-632291、
  SnS mp-554134、層状 BN mp-685145、H8 mp-730101 など。試すこと: ES の基底に p、すき間に ES を足す、H の基底（EH2 の s,p）。
  層状の結晶にも ES のしきい値 2.0 a.u. を使う（結合のつながりの次元で分子・層状を判定する）
- **【未着手】`wcsmear = false` と tf32 で `hgw` が Σc の途中で落ちる**（2026-10-02、MgO、メッセージ無し）。`wcsmear = false` は既定ではないが、比べるときに使う

## 試験 {#tests}

- **【未着手】板の ES の試験をサンプルにする**（2026-10-08）。PbS の板（mp-561320）で、ES の置き方により MLO が通る・通らないの差がはっきり出た。
  QSGW80 は 5 時間かかるので、LDA か小さい板で同じ差が出る形を探す

## 文書（ecaljdoc） {#docs}

- **【判断待ち】残した朱の印**: [MLO](./mlo) の 2 か所（損失の履歴の「0.25 Ry」、「12 系」の注）。コードやデータで確かめきれていない
- **【未着手】数式の本文からのはみ出し**（2026-11-14）。幅 1300px で kBT・pmt・mlo・gw、960px ではさらに basic・basic2・mlo_gwsc。`aligned` で折り返す
- **【未着手】メニューの見出し「GPU version of QSGW」を日本語に**
- **【未着手】同じ説明の二重化の残り**（2026-11-14）。メインドキュメント（README_tutorial）とほかのページ（`save` の文字、lmchk、SOC の手順など）、kBT・mlo の中の重複、
  Samples の README の古い記述
