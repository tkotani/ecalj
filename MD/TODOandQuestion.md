# TODO と質問、やったこと

直すべき点は見つけてもその場では直さず、ここに書く（user 2026-10-01）。メンテナに決めてほしいことも、やったことも、ここに書く。
入口は [CLAUDE.md](../CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）。片付けたものの中のノウハウは [past_log.md](past_log.md)。
細かい経緯と数値は研究ログ [research_log.md](research_log.md)。

書き方: 各項目に見つけた日を付ける。済んだら「やったこと」へ移し、済んだ日とコミットを書く。

---

## 1. TODO（直すべき点）

### 進める順序（2026-10-02 05:24 に決めた。user「TODO の手順はよく考えて順序立てて」「まかせる」）

依存の向き: 対称性の整理（S1〜S6）が、MLO の最大局在化（Sakuma 型の拘束に操作が要る）と AF の項目の土台。最大局在化の結果で、マグノンの窓と
Wannier とのずれの項目の判断が変わる。だから対称性 → 最大局在化 → マグノン・ずれの順。長い試験の間に、独立した小さな項目を挟む。

| 順 | 項目 | なぜこの位置か | 確かめ | 状況（2026-10-02 06:57） |
| --- | --- | --- | --- | --- |
| 1 | 対称性 S1（分ける、等価変換） | 以後の対称性の段の土台 | 172 入力の lmchk が S0 と同じ、試験の組 | **済み** `41c0e3cc7` |
| 2 | 小さい独立の直し: `job_mlo_soc` の空の spin2、`m_tetrakbt` の使われないルーチン、`auto_creplot.py` | 1・3 の試験の待ち時間に | mlo、kBT の組 | 前二つ**済み**。`auto_creplot.py` は移すか捨てるか要判断 |
| 3 | 対称性 S2（GW 側の `mptauof` の重複を外す） | S3 の前に使う側を一本に | gwall がビット単位で同じ | **済み**（最小の形）`e5d08f008` |
| 4 | GPU の build の module の循環 | kt1・kr7 で回せる、手元と並行 | kt1・kr7 の clean build と試験 | **済み** `a7f64752e`・`f553b5222` |
| 5 | 対称性 S3（`symmetry.json` を読む口）→ S4（純粋な並進）→ S5（AF、`AFsymmetry/NiO` の pwmode も）→ S6（既定に） | 設計どおり一段ずつ | 各段の表（`MD/symmetry_spglib.md` §4.7） | S3〜S6 **済み**（S6: lmf が同梱の spglib で求める、user の判断 2026-10-02。2026-10-02 08:34 から 3 台で試験） |
| 6 | MLO の最大局在化（Python で試作、Fe・Ni） | 5 の操作を使う | Ω、U、マグノン | 試作**済み**（`mlo_maxloc.py --sym`）。Ω まで。U・マグノンは未 |
| 7 | マグノンの既定の窓、MLO と Wannier のずれ | 6 の結果で判断 | Fe・Ni・FeCo | 3 物質で比べた（07:10）。Δ = 6 eV を勧め、§2 で判断待ち |
| 8 | MLO の模型の残り（EH2 の崩れ、§9 の目安の値の測り直し、空格子球の自動化） | 計算機で裏で回せる | MATERIALS | 原因の確かめと測り直し、`ctrlg_addes.py` の直し**済み**。EH2 の直し方は要判断 |
| 後 | `sugw` のメモリ、`hgw` の残り、MP の API と GW1500 の選定 | 大きい、または外の事情 | — | 未 |


### コード

- **MLO: EH2 s,p を単体の金属（Cu・Ni）に足すと、Γ–X の 2〜3 点の k だけで模型が崩れる**（2026-10-01、`~/work/mlocheck_eh2cat/{Cu,Ni}`、ページ https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH の図 4）: その k で MLO の帯が E_F + 0.3 eV に集まり、DFT の帯が抜ける（Cu 0.012 → 0.802 eV）。同じ原子の EH と EH2 がほぼ一次従属になって `Hreduction` の重なりが特異に近い、と疑っている（未確認: その k の重なり行列の固有値を見る）。gwinit の `!` 付きの `mlo_lm2` の行は Cu・Ni には書かれない（既定は基準 1 で EH2 なし）ので、既定には影響しない。EuO（Eu に EH2 s,p）は模型を作る所で止まる: `m_hreduction` の NormalizationCheck でシードのノルムの減りが 1 % を超えた（band 26、−1.3 %。EH と EH2 の両方をシードにすると `zhev_tk4` が一次従属に近い向きを落とすため）。Cu・Ni も同じ原因で、減りが 1 % 未満なので規格化し直して進み、特定の k で崩れているのかもしれない（未確認）
  2026-10-02 06:22 確かめた（`~/work/eh2cu`、前日の DFT から今の版で `job_mlo`）: 帯の道筋の x = 0.429 で MLO の重なり O(k) の最小固有値が **負**（−1.26×10⁻⁴、
  道筋の中央値 3.2×10⁻⁴）、その近く x = 0.381 で帯が最大 8.2 eV 外れる。実空間で打ち切った O(k) のフーリエ和が、EH と EH2 がほぼ一次従属だと正定値で
  なくなり、`m_mlo_ham` の `calc_ham_eigen` の `zhgv`（Cholesky）が壊れる。案: (i) 対角化で O(k) の小さい固有値の向きを捨てる（正準直交化。帯の数が k で
  変わるので出力の形を決める要あり）、(ii) 模型を作るとき（`Hreduction`）一次従属に近い EH2 を自動で外す。既定（基準 1、EH2 なし）には影響しない。どちらにするか要判断
- **`sugw`（`lmf --jobgw=1`）のメモリ**（2026-10-01）: `GEIGpart` が IPW の重なり行列 `ppovl(ngp,ngp)` と LU 用の写し `ppovlLU` を
  各ランクで持つ。1 ランクあたり 32·ngp² バイト、全体は並列数に比例。ngp ≈ V·Q³/6π²（Q = `QpGcut_psi`）なので胞の体積の 2 乗で増える
  （Rb8 の 4410 Å³ で 1 ランク 35 GB）。案: (a) 投入前に ngp から並列数を決める、(b) O は正定値なので Cholesky の因子だけ持つ（半分）、
  (c) O·x を FFT で作り反復法で解く（ngp に比例）
- **`auto_mpquery.py` が今の Materials Project で動かない**（2026-10-01）: MP の API が新しい ID（`mp-aaacpdie` の形、古い番号の 26 進）を返し、
  手元の `mp_api` の検証（pydantic）が止まる。`mp_api` を上げるか、REST を直接読む（requests、`x-api-key`、User-Agent が要る。urllib の既定は 403）
- **GW1500 の選定に構造の確かめを入れる**（2026-10-01）: 別の化合物から副格子を抜き出した MP の項目（ICSD の備考が「〜 part」、
  1 原子の体積が同じ組成の 2 倍以上など）が PBE のギャップ > 0 で選ばれていた（`ecalj_auto/GW1500_status.md` 表 6 の `INVALID_STRUCTURE`）。
  `auto_mpquery.py` で弾くか、印を付ける
- **AFTEST（`mmtarget.aftest`、例は `Samples/AFfixMMOM`）の残り**（2026-10-01、研究ログ 2026-10-01 朝 06:46 と表 06:46-1）: 場をスピン 1 にしか入れていなかった誤りは直した
  （`6eff2df53`）。残り: (a) モーメントを下げる向きで更新が行き過ぎる（利得 2 固定と 0.1(m − m_t)² の項。NiO 1.28 → 1.0 は 70 反復で、途中で符号が反転）。
  応答を割線で見積もる更新にするか。(b) 対がサイト 1・2、ブロック 1・2 の決め打ち（NiSe だけ 6 ブロック）。`AF=` の印や種の名前から対を決める。
  (c) `m_ldau_init` が lmf の起動のたびに場を一歩更新して `mmagfield.aftest` を書き直す（`job_band` でも）。(d) ecaljdoc の `UsageDetailed.md` の「直す必要がある」を、
  使い方（LDA+U のブロックが要る、U = 0 でよい、`SYMGRPAF`、目標のモーメントの定義）に書き直す（2026-10-02 済み: ecaljdoc の UsageDetailed.md に移し、`Samples/AFfixMMOM/README.md` は参照だけに）

- **MLO と Wannier のずれの追跡**（2026-10-02、`MD/wannier_vs_mlo.md`）: Ni d の cRPA の U が MLO で 25 % 低い（d の部分空間の切り出し方の違い）、Fe のマグノンの q ≥ 0.4 で MLO が 1.5〜2.5 倍高い。実験（Fe のマグノン分散）との比較と、MLO の窓（`mlo_delta`、`mlo_w`）への依存を測る
- **対称性: Fortran の中で「見つける・導く・使う」を分け、spglib を入れる**（2026-10-02、`MD/symmetry_spglib.md` §4.7 の S0〜S6）。S0 済み（`f7c382bed`、食い違い 0）。超格子は暫定で閉じた部分群（Si の慣用胞で 24）、純粋な並進を使うのは S4。2×1×1 の Si の GW は遅いだけだった（§5）
- **MLO のマグノンの既定の窓**（2026-10-02、研究ログ 04:17）: 窓を広げる（Δ = 6 eV）と d の MLO が局在し、η ≈ 1、高い q のマグノンが Wannier 版に近づく。`job_mlo_magnon` 用に既定を変えるか、`Fe_mlo_magnon` の入力を変えるか（参照を作り直す）。Ni、FeCo でも確かめる
- **MLO の部分空間の中での最大局在化**（user 2026-10-02 未明）: MLO を作った後に Marzari–Vanderbilt で局在させる（バンドと cRPA の p_kn は変わらない）。MLO は直交していないので、まず Löwdin で直交化し、そこから最小化する。対称性は Sakuma（PRB 87, 235109）のように U(gk) = D(g) U(k) d(g)⁻¹ で拘束する（MLO は MTO と同じに回る、`rotmatMTO`）。`mlo_spread.py` の M(k,b) が材料。比べる量: Fe・Ni の d で Ω（生、Löwdin のみ、Löwdin + MV）、マグノン、U
  2026-10-02 06:20 試作 `SRC/exec/mlo_maxloc.py`（拘束なし）: Ni d は Ω の 97 % が Ω_I で MV は効かない。Fe spd では拘束なしの MV が対称性を破る（研究ログ 06:20）。06:31 に Sakuma の拘束を入れた（`--sym`、1 原子・symmorphic だけ）: Fe t2g 7 % 縮むだけ。残り: 多原子・非 symmorphic、局在させた基底での cRPA の U とマグノン
- **`hgw` の残り**（2026-04 の統合の残タスク、2026-10-01 に past_log.md §3.2 から）: ノード内の W の共有を `MPI_Win_allocate_shared` で（メモリの重複を減らす）。
  `hsfp0_sc` の Sx（`--job=1`）・core の交換（`--job=3`）も `hgw` に入れる（時間は小さいので優先度は低い）


- **MLO の模型の残り**（2026-10-01。既定は基準 1・2・3 の形に決まった: ecaljdoc mlo §1・§9、`MD/handover.md` §5）:
  (a) 基準 3（空格子球）の自動化。案は空隙の半径（最も近い MT 球の表面までの距離）が 3 a.u. を超える所に置く `SRC/exec/ctrlg_addes.py`（未コミット・未検証）。
  Bi₂Te₃ の空隙の半径は 2.61 a.u.（中心 (½,½,½)、最大の点で 2.68）で 3.0 より小さい。前の 3.56 は誤り（2026-10-02 06:26 に確かめ、道具の像の数を面間隔から決めるように直した）。
  SiO₂ は手で置いて 0.001 eV（研究ログ 2026-10-01）。
  (e) 検査（`mlo_bandcheck.py` の CHECK）で、基準 1 の Sn・AlSb・Bi₂Te₃・InSb・Cu が 1 点だけ 0.1 eV を少し超えて FAIL（rms は 0.01〜0.03）。
  線形独立性の崩れではない（重なりの最小固有値は 3e-3 程度）。帯の交差の付近の形が違う。直すか、目安を見直すか
- **`Samples/AFsymmetry/NiO` は `pwmode = 1` で `symgrpaf`**（2026-10-01）: k 点を 4³ にすると `rotwave: q+G rotation error (We have to set PWmode=11 for symgrpAF)`
  で止まる。試験の 3³ では通っているだけ。`pwmode = 11` にして参照を作り直すか

### 試験と入力

- `MLOsamples/RuO2` の保存してある `rst`・`dmats` は、`pwmode = 11` の LDA+U の誤り（2026-03-30〜09-30、`2498e5283` で修正）の時期に作ったもの。
  作り直すか（GdCo5・SmP は 2026-09-30 に作り直した）
- 試験の入力の温度（`t_tetrakbt = 262`、`t_sigmaw = 0`）をテンプレート（300/300）に揃えるか。揃えると gas_gwsc・fe_gwsc などの参照が動く
- `SRC/exec/auto_creplot.py` がまだ旧形式の `ctrl.<sname>` を書き換える（2026-10-01）。ctrlg に直すか、使われていなければ trash へ。（2026-10-02: `auto_job_mp.py` が `import creplot` で使い、`ecalj_auto/auto/creplot.py` に写しがある。GW1500 の回し直しは `ctrlgenToml.py` と `gwscconv` を直接使っていて通らない）

## 2. 質問（メンテナに決めてほしいこと）

- **MLO のマグノンの窓を Δ = 6 eV にするか**（2026-10-02 07:10、研究ログ 07:10）: Fe・FeCo・Ni で、Δ = 6 eV（w = 2）は Wannier のマグノンから 0.02 eV 程度、
  既定（2, 2）は Ni では良いが Fe・FeCo の高い q で 2 倍近く高い。`job_mlo_magnon` が窓を Δ = 6 eV にするか、試料の入力に書くか。`Samples/Magnon/Fe_mlo_magnon` の参照の作り直しを伴う
- **ecaljdoc の古い文書**（`BackUp/`、`ecaljdetails/` の LaTeX・PS、2019 年以前）は、trash に移した `TOOLS/checkmodule`・`TOOLS/ModuleCodingSample` などを参照している。ecaljdoc の側も同じ決まり（trash へ、要点は過去ログへ）で片付けるか（2026-10-01）
- **GW1500 の `INVALID_STRUCTURE`（12）・`SUSPECT_STRUCTURE`（4）を集合から外すか**。いまは注記だけ
- **ビルドの生成物が入ったコミット `a0c7a7300`（78 MB）を、push の前に履歴から消すか**。消すと以後 674 コミットのハッシュが変わり、
  研究ログなどに書いたハッシュが合わなくなる。消さなければ、公開のリポジトリの履歴に 78 MB が残る
- **Materials Project の API キーの作り直し**（user の作業）: 2025-04〜2026-09 の公開リポジトリの履歴に入っていた
- **trash を空にするか**: t14 の `~/ecalj/trash`（旧ツリー 58 GB ほか）、kt1 の `~/ecalj/trash`（255 GB）。kt1 は空にするまで `/` の空きが 49 GB のまま
- 試験に使った古いツリー（mic の `~/ecalj_test0928`、kt1 の `/mnt/data1/ecalj_test0930b`・`0930c`）も trash に入れるか
- ブランチ `fix-idu10`（main にマージ済み）を消すか
- push: dev・rel とも、t14 の main より 674 コミット遅れ（2026-10-01）

## 3. 実行中

- kt1 `/mnt/data1/gw1500_rerun/run3`: GW1500 の NOTCONV の残り（fp32、5 本）。2026-10-01 04:40 に約 17 時間の見積もり

---

## 4. やったこと（新しい順）

### 2026-10-02

- 対称性 S6（2026-10-02 08:34、user「lmf でつくればいい」）: spglib 2.6.0 の C を同梱し、lmf・lmchk が操作を求めて `symmetry.<sname>.json` を書く（`89654cc96`）。172 入力で Python 版と同じ
- ecaljdoc mlo §9 の式 (12) の目安の値を測り直した（06:55、`~/work/ovlp_20261002`、63 物質の MLO の段だけ）: 規格化の後は最小 0.20〜0.42、中央値 0.24〜0.46。FAIL の 10 物質は前と同じ（ecaljdoc `mlo.md`、TODO から外した）
- 対称性（`MD/symmetry_spglib.md` §4.7）: S0 spglib との照合（`f7c382bed`、172 入力で食い違い 0）、S1 分割（`41c0e3cc7`）、S2（`e5d08f008`）、
  S3 `symmetry.<sname>.json` を読む口（`91f9bdcaf`）、S4a `mptauof` が渡された並進を使う（`bfaae451b`）。S4b（純粋な並進を操作に）は試験中
- GPU の build の module の循環を切った（`a7f64752e`）、CMake の回避策を外した（`f553b5222`）。kt1・kr7 でまっさらなビルドが通った（TODO から外した）
- `m_tetrakbt` の使われないルーチン（`8217b212b`）、SOC の MLO が空の spin2 を書く件（`0a9f1b092`）（TODO から外した）
- MLO を実空間で規格化（`0e33fbd87`）、`mlo_spread.py` と比較の記録（`ec4361edc`）、比較の一式とタグ `last-wannier`（`f1de3817a`）、cRPA の試験を MLO 版に（`e22bb9075`）、Wannier・AHC・lmfham2 を外した（`6d8b9f05b`、`dfffe0dfd`）。`MD/wannier_vs_mlo.md`、ecaljdoc mlo §6
- MLOsamples などの古い作業ファイル 68 本（`ctrlp.*`、`lmfham2parameters.check`、`out_lmfham1`、`bandplot_MPO.*` など）を trash へ（past_log 表 1）
- `TOOLS/sync_ecalj_src.sh`: 送り先の `SRC/subroutines`・`main`・`exec` にあって HEAD に無いファイルを `trash/` へ移す（Wannier を外したとき kt1・kr7 に 30 本残って CMake が拾った）

### 2026-10-01

- Samples の旧形式の入力を trash へ（user の指示）: TestInstall の `ctrl.*`（`5b0c3207d`）、残りの Samples の `ctrl.*` 28 本（`f2f79a380`）、MATERIALS の `GWinput`・`ctrlgenM1.ctrl.batio3` と `MLOsamples/Al2O3_Cr/CASE1ok`〜`CASE5ok`（`b9cf99a4a`）。中身の要点は past_log.md 表 1
- MLO の模型の検査（`mlo_bandcheck.py` の CHECK、`job_mlo`・`job_mlo_soc` が最後に回す、`mlo` が `MLO_ovlpmin.dat` を書く）、
  窓の基準の E_F を `--efermi=` のファイルから取る直し、`Samples/MATERIALS` の 63 の ctrlg の `[mlo]` を今の gwinit で書き直し（研究ログ 2026-10-01 21 時）。
  試験: mlo 45、MLO-QSGW 5、install 64、inputs 176 が t14 で PASSED
- MLO の模型を基準 1・2・3 に整理した（user と決めた。ecaljdoc mlo §1・§4・§9、`MD/handover.md` §5、研究ログ 2026-10-01 16〜20 時）:
  半内殻の局所軌道は帯の上端で自動（E_F − 8 eV より上は EH と入れ替え、−17〜−8 eV は加える、`d759110e4`）、gwinit は f と `!` 付きの `mlo_lm2`（基準 2）を書く、
  MLO-QSGW の凍結した模型の本数の誤り（`2053bb5ea`）、`mlo_bandcheck.py` の窓を [VBM − 8, CBM + mlo_delta] に（`bd1b893ec`）、
  `job_mlo_soc` がスピン軌道ありの DFT のバンドも描く（GaAsSoc の −0.1 eV は比べ方の誤り、`dc66439c5`）。Samples/MATERIALS の結果のページ
  https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH（`Samples/MATERIALS/mlocheck/gallery.py`）。試験: inputs 176、mlo 25 試料、MLO-QSGW、Fe_mlo_magnon が t14 で PASSED
- AFTEST の修正（`aftest-fix`）を main にマージ（`6eff2df53`、user「マージできるよね」→「良い方を壊さないように」）。AF の対称性ありの計算は ehf・sev・モーメント・場がマージ前と同じで、ehk だけが制約の下の全エネルギーになった。afsym 4・install 64 が PASS。例を `Samples/AFfixMMOM`（user の命名、`NiO_afsym`・`NiO_noafsym`、試験の組 `affix` 12 件 PASS）に置いた（`ca1e4c7fa`）
- ecaljdoc の TOML の流れの最短の手順は、`manual/README_tutorial.md` の GetStarted（Step 0〜6: POSCAR → ctrls → `ctrlgenToml.py` → `lmfa`・`lmf` → バンド → `gwsc`、Step 2-Migration、`--ctrlg:` の上書き）に既にあった。TODO から外した。ecaljdoc `manual/mlo.md` に既定の模型が外れる 2 つの型と表 M1 を書いた（`422dc92`、未 push）
- `Samples/MATERIALS` の 65 物質で LDA と MLO の自動の模型を回した（InAs/GaSb n10 の 40 原子は user の判断で回さない）。結果の表は `Samples/MATERIALS/MLOcheck_20261001*.tsv`、一覧のページ https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH、研究ログ 2026-10-01 朝 07:16
- `job_mlo`: ctrlg の `[ham] so = 1` なら `job_mlo_soc` を案内して止まる（`mlo` がハミルトニアンの NaN で止まっていた。GaAs_so）。`--ctrlg:ham.so=` の上書きは尊重。`mlo_bandplot.py` は空の spin2 を描かない
- `--cls`: `m_clsmode_finalize` に `ndimh`（`m_igv2x` の最後の k 点の値）でなく `nbandmx` を渡す（`vcdmel` の重みの並びと `dostet` の読み方を揃える）。CrN（`Samples/TestInstall/crn`）で APW なし・pwemax 2 と 5・4×4×4（ndimh が 182〜188 と変わる）のどれも `dos-vcdmel.crn` が一致（ずれるのは DOS の窓より上の帯だけだった）。ブランチ `cls-nbandmx` をマージ（`7ca6f1caa`）。試験は `~/work/clstest_20261001`
- `InstallAll.py`: build の後に、bindir の中で「この ecalj の木を指していて行き先の無いリンク」だけを消す（`remove_dangling_links`）。`SRC/exec` の行き先の無いリンク（消した `SRC/exec/build/` を指す 12 本）は張らず、trash に移した。これが `~/bin/mlo` などを一時的に死んだリンクで上書きしていた（CMake の deliver が build の後に戻していた）。SRC/exec のエディタの一時ファイル（`~` で終わる、`#`・`.#` で始まる）はリンクしない。t14 の `~/bin` の 58 本は次のインストールで消える（一時の bindir で試験）
- ecaljdoc `manual/spectrum.md` に、メッシュの外の q（`QforGW`）の注意 3 点（`EMAXforGW` が必須、窓を変えたら交換から、`epsWVR` の行）を書いた（past_log.md §5 から、ecaljdoc の未 push のコミット）
- AFTEST を調べた（研究ログ 2026-10-01 朝 06:46）: afsym ではモーメントを目標に保てるが、表示の ehk に −uhx·m_d、ehf に −2·uhx·m_d が残る。afsym なしでは誤り（サイトの電荷が分かれる）。修正はブランチ `aftest-fix`
- GW1500: kr7 で fp32 と TF32 を同じバイナリ・同じ入力で比べた（8 物質、02:03〜06:24）。最終のギャップの差は 1 meV 未満。5 月との 0.29〜1.48 eV の差は精度ではなく、5 月の振動と設定の違い（`GW1500_status.md` §5.1 の表 7）。GOOD 1120 を精度の理由で見直す必要は無い
- `SRC/subroutines/m_tetrakbt_BUGREPORT.md`（2026-06、直し済みの不具合の報告）の要点を past_log.md §13 に移して trash へ。空の道標 `SRC/TestInstall_is_moved_to_under_ecaljSamples` も trash へ
- `.gitignore` に `/build/`（VSCode の CMake 拡張が最上位に作る）
- `MD/module_map.md`（生成物）と `TOOLS/module_map.py`: 主プログラム → 入口の module、module の階層（226 module、最大 29 段、`m_lmf` が頂点）、
  依存と被依存の数、module の冒頭のコメント。`ecaljclaude.md` から参照。GPU の build の module の循環が 1 つ見つかった（上の TODO）
- `GetSyml/README.md`・`StructureTool/README.md`・`Samples/EPS/EPS_GaAs/README_eps.md` を今の形に書き直した（入力は `ctrlg.<sname>.toml`、
  `getsyml --nobzview`、StructureTool の各スクリプトの向きと出力のファイル名、EPS は `[gw] QforEPS`・`QforEPSau`・`n1n2n3`・`[product_basis] pb_lcutmx`、
  EPS の出力の列）。`getsyml` の使い方の表示も（`-nobzview` → `--nobzview`、`ctrl.nio` → `ctrlg.nio.toml`）
- `ctrlgenToml.py`: nspin=1 のとき原子表の IDU/UH/JH をコメントにして書く（LDA+U は nspin=2 が要る）。`Database/Ce`（nspin=1、idu=12）が
  09-30 の idu の修正（`78475ea3d`）から lmf で止まっていた（`LDA+U must be spin-polarized!`）ので同じく直した（`2cbeb93d5`）
- `Samples/MATERIALS` を仕分けた: 構造のデータベース（62 物質）を `Database/` に展開（`ctrls` と、GW・MLO の節つきの `ctrlg`、`lmchk` で全部読める）。ほかの旧形式の 30 項目は trash（624 ファイル）。残したのは `La2CuO4`・`InAsGaSb`・`BaTiO3`。拾ったノウハウは past_log.md §9
- `README.txt` と `.org` を Markdown に（user「README.txt とあるのは md 形式に。org もそう」）: `StructureTool/README.txt`、`Samples/EPS/EPS_GaAs/README_eps.org`、`Samples/MATERIALS/{LaGaO3_relax,MLOsamples,NiSe_aftest}/README.org`、`Samples/MATERIALS/Si_doping_sample/Memo_bgcharge.org`（`git mv`）。pandoc は `_` を下付きに、`--` をダッシュに変えるので使わず、原文を保つ変換（見出し・`#+TITLE`・`#+begin_src`・コマンド行と設定の断片をコードブロックに）
- `GetSyml/`・`StructureTool/` の古い例（旧形式の `ctrl.*` 86 本、`syml.*`、鉱物の POSCAR 135 本）と使わないスクリプトを trash へ（257 ファイル）。本体（`getsyml`、`vasp2ctrl`・`ctrl2vasp`、`viewvesta`、`refineposcar.py`、`superlattice/`）は残し、動くことを確かめた。past_log.md §12
- `Doxygen/` を trash へ（user「そうしよう」）。コメントを Doxygen 形式に揃えることはしない。作り直し方は past_log.md §11
- `TOOLS/` の古い道具（約 60 項目、632 ファイル）を trash へ。残したのは `samples_tests.sh`・`sync_ecalj_src.sh`・`ozbench/`。中身は past_log.md §10。`diffnum` は試験で今も使うが、使うのは `SRC/exec/pylib/diffnum0.py`（TOOLS の版は古い）
- Claude の個人メモリから、引き継ぐ価値があり今も正しいものを [handover.md](handover.md) に写した（user「メモリの内容はパッケージに入らないので MD/ に。混乱を招くものは良くない」）。
  写す前に 5 点をコードと照らした（rel の既定ブランチは `main`、`master` は無い、など）
- `README.md` を `MD/README.md` へ（user「Claude に読ませて、人間は Claude から情報を取る構造にする。README を人間に読ませるのは好ましくない」）。
  最上位の `README.md` は数行の案内だけ（GitHub の表紙が空にならないように）。`ecaljclaude.md` の方針の文言を直した
- 開発の文書を `MD/` に（`ecaljclaude.md`、`TODOandQuestion.md`、`past_log.md`、研究ログ `research_log.md`、`bf8985687`）
- ecalj の片付けの 3 回目: `PHASE1B_REFACTOR.md`、`HIGHLIGHTS_2026-06_09.md`、`FiniteT_and_QPE_HOWTO.md`、`ecaljdoc_drafts/`、`jobauto/`、
  `SRC/exec_legacy/` を trash へ。中身は past_log.md §3.2・§4.3・§5〜§8 に整理し、README と ecaljdoc（ForDevelopers・kBT・README_tutorial）の参照を直した。
  最上位の `MATERIALS/` は trash ではなく `Samples/MATERIALS/` の下へ移した（user「いったん Samples の下へ」、624 ファイル、`git mv`）
- ecalj の片付け（user「いらないものは trash へ。ノウハウは過去ログへ。重複は整理してから trash へ」）:
  ビルドの生成物 1853 ファイル（`3f0771f2f`）、最上位の打ち込み用と古いスクリプト、`.refactor_notes/`、`SRC/exec/BK`、
  MLOsamples の古い試行、`TestInstall/TESTunused`、`*.bk` などを trash へ。ノウハウは [past_log.md](past_log.md)
- 最上位の `SRC/BK`・`SRC/execgfortran`・`SRC/execAHC`（181 MB）と、GW1500 の古いスクリプトを追跡から外して trash へ（`116254e1d`）
- kt1 の GW1500 の作業ファイルと退避（255 GB）を `~/ecalj/trash` へ（`4f9332d98`）
- ecalj・ecaljdoc の控えのクローンを外付けドライブ `/media/takao/TAKAOMINI`（`ecalj_20261001`・`ecaljdoc_20261001`）に
- Materials Project の API キーを追跡から外した。`<ecalj>/MaterialProject.key`（無視される）と見本 `.example`、`pylib/mpkey.py`（`290397b34`）
- GW1500: 物質ごとの注記の表 `ecalj_auto/gw1500_notes_20261001.tsv`、分類と条件の違い（`GW1500_status.md` §5、`40ffd89e8`）。
  5 月に収束した GOOD はやり直さない（run3 のキューを NOTCONV だけに）。Rb8（mp-1179832）は不正な構造と分かり止めた（`e616f55fb`）
- `gwscconv`: ギャップの無い反復（金属）は固有値の変化で収束を判定（`ecfac2c6b`）
- kt1・kr7（nvfortran GPU）・mic（ifx）で試験の組 inputs・mlo・afsym・install・samples がすべて PASS（HEAD `74ba72dad`）

### 2026-09-30

- `fix-idu10` をマージ（`78475ea3d`）: `idu = 10 + mode` が `sigm` の無いとき二重計数なしの LDA+U になっていた（2023-09 から）。
  MLOsamples の SmP・GdCo5 の LDA+U を作り直した（`e547b81e0`、`a7b19be75`）
- GW1500: FAILED の 147 物質を fp32 で回し直し、144 が収束（`d00b9778e`）
- Samples/Legacy を片付けて組み直し（AtomDimer/N2 など）、ecalj と ecaljdoc を clone し直した
