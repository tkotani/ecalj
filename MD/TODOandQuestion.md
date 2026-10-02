# TODO と質問、やったこと

直すべき点は見つけてもその場では直さず、ここに書く（user 2026-10-01）。メンテナに決めてほしいことも、やったことも、ここに書く。
入口は [CLAUDE.md](../CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）。片付けたものの中のノウハウは [past_log.md](past_log.md)。
細かい経緯と数値は研究ログ [research_log.md](research_log.md)。

書き方: §1 はまだのものだけ（各項目に見つけた日）。済んだら §2 へ移し、済んだ日とコミットを書く（user 2026-10-02「まだのものを冒頭に、終わったものは 2. TODO（済）に」）。
計算機ごとの片付け（trash の削除、ディスクの空き）はここに書かない（user 2026-10-02「ローカル情報はコミットに入れなくていい」）。

---

## 1. TODO（未）（2026-10-02 14:35 に組み直した）

### MLO

- **Löwdin を標準にした後のメッシュの外の悪化**（2026-10-02、研究ログ 12:16 の表 12:16-1）: Fe の 4s 帯の底（−8〜−3 eV で 0.03 → 0.096 eV）、E_F 近くの 2H-SiC（0.035 → 0.051）・
  MnO（0.007 → 0.030）。MLOsamples の Fe 系の検査も 0.217 → 0.284。原因を調べる（その k・帯、H̃(R) の切り方）
- **`mlo_bandcheck.py` の判定 (3)**（2026-10-02）: 模型の O が 1 になり、重なりの最小固有値（`MLO_ovlpmin.dat`）は常に 1。`mlo` の出力の
  "Smallest eigenvalue of the normalized raw overlap"（メッシュ上の直交化の前）を判定に使うように替えるか
- **EH2 を足した模型の崩れ**（2026-10-01、`~/work/eh2cu`、Cu・Ni、EuO）: 生の模型では実空間で打ち切った O(k) が正定値でなくなり `zhgv` が壊れた。
  Löwdin の模型は O(k) を内挿しないので様子が変わるはず。測り直してから直し方（(i) 正準直交化、(ii) 一次従属に近い EH2 を自動で外す）を決める（要判断）
- **η ≠ 1 の原因**（2026-10-02）: Goldstone の条件の倍率が Löwdin で Fe 1.24、FeCo 1.28、Ni 1.75。Ni が大きい理由は未確認
- **MLO と Wannier のずれ**（2026-10-02、`MD/wannier_vs_mlo.md`）: Ni d の cRPA の U が Wannier より低い（Löwdin で 2.90 対 3.78 eV、部分空間の違い）。
  Fe のマグノンは Löwdin で近づいた（q ≤ 0.3 で 2 割高いのは残る）。実験（Fe のマグノン分散）との比較
- **最大局在化（MV）の試作の残り**（2026-10-02、メインにはしない）: `mlo_maxloc.py --sym` は 1 原子・symmorphic だけ。多原子・非 symmorphic、局在させた基底での
  U とマグノン。`mlo_cmlo_transform.py lowdin` は恒等変換になった。続けないなら試作の道具をまとめて trash へ
- **空格子球の自動化（基準 3）**（2026-10-01）: `SRC/exec/ctrlg_addes.py`（`8c158ee96`、像の数を面間隔から決める直しは済み）を ctrlg の生成に組み込むか。
  Bi₂Te₃ の空隙は 2.61 a.u. で 3.0 未満、SiO₂ は手で置いて 0.001 eV
- **検査の FAIL の残り**（2026-10-01）: 基準 1 で Sn・AlSb・InSb が 1 点だけ 0.1 eV を少し超える（Löwdin の模型で Bi₂Te₃・Cu は通るようになった。表 M8）。
  AlN・MgS・MgSe・MgTe・SiO₂ は誤差が大きい。直すか、目安を見直すか
- **`Samples/AFsymmetry/NiO` は `pwmode = 1` で `symgrpaf`**（2026-10-01）: 4³ にすると `rotwave: q+G rotation error` で止まる。`pwmode = 11` にして参照を作り直すか

### 有限温度

- **有限温度の四面体法（`t_tetrakbt > 0`、`m_tetrakbt`）と従来の T = 0 の四面体法（`t_tetrakbt = 0`、`tetwt5`）の関係を確かめる**（2026-10-02 16:02、user）:
  T → 0 で有限温度の重みが従来の重みに一致するか（同じ k メッシュ・同じ対称性で、χ0 の虚部と QP のエネルギーを比べる。T を 300、100、30、10 K と下げる）。
  Σ 側の Fermi–Dirac の幅（`t_sigmaw`）と χ0 側の温度の組み合わせも。きっかけ: heavy の `nio_gwsc444`（`t_tetrakbt = 0`、`t_sigmaw = 262`）の QPU が
  2026-10-02 の版で 15 meV ずれた（kr7・kt1 で同じ値）。原因は確かめ中（従来の対称性の探し方で回し直して切り分け）

### コード

- **HEAD に残っている生成物らしいもの**（2026-10-02、旧 `a0c7a7300` で入り今も追跡）: `SRC/.#Memo4rotation`、`SRC/exec/cmake_install.cmake`、`hello.py`、`platform`（800 KB）、
  `lmf2.py`・`lmchk.py`・`pylmfa`・`pysample`・`ohtaka`・`epsPPd`・`epsPPsaito`・`job_senefbz`・`readeps_dig2.py`・`auto_kauto.py`。使われているかを確かめて trash へ
- **`sugw`（`lmf --jobgw=1`）のメモリ**（2026-10-01）: `GEIGpart` の `ppovl(ngp,ngp)` と `ppovlLU` を各ランクで持つ（32·ngp² バイト、胞の体積の 2 乗）。
  案: (a) ngp から並列数を決める、(b) Cholesky の因子だけ持つ、(c) O·x を FFT で作り反復法で解く
- **`hgw` の残り**（2026-04、past_log.md §3.2）: ノード内の W の共有（`MPI_Win_allocate_shared`）、`hsfp0_sc` の Sx・core の交換も `hgw` に（優先度は低い）
- **AFTEST の残り**（2026-10-01、研究ログ 2026-10-01 06:46）: (a) モーメントを下げる向きで更新が行き過ぎる（割線で見積もるか）、(b) 対がサイト 1・2、ブロック 1・2 の決め打ち、
  (c) `m_ldau_init` が lmf の起動のたびに場を更新して `mmagfield.aftest` を書き直す（`job_band` でも）
- **`SRC/exec/auto_creplot.py`**（2026-10-01）: 旧形式の `ctrl.<sname>` を書き換える。`auto_job_mp.py` が使う。ctrlg に直すか、使わないなら trash へ

### GW1500 と Materials Project

- **`auto_mpquery.py` が今の MP で動かない**（2026-10-01）: 新しい ID の形で `mp_api` の検証が止まる。`mp_api` を上げるか、REST を直接読む
- **GW1500 の選定に構造の確かめを入れる**（2026-10-01）: 副格子を抜き出した MP の項目が選ばれていた（`ecalj_auto/GW1500_status.md` §5.3 の表 8・9）。`auto_mpquery.py` で弾くか印を付ける

### 試験と入力

- `MLOsamples/RuO2` の `rst`・`dmats` は `pwmode = 11` の LDA+U の誤りの時期（2026-03-30〜09-30）に作ったもの。作り直すか
- 試験の入力の温度（`t_tetrakbt = 262`、`t_sigmaw = 0`）をテンプレート（300/300）に揃えるか。揃えると gas_gwsc・fe_gwsc などの参照が動く

### メンテナに決めてほしいこと

- **GW1500 の `INVALID_STRUCTURE`（12）・`SUSPECT_STRUCTURE`（4）を集合から外すか**。いまは注記だけ
- 試験に使った古いツリー（mic の `~/ecalj_test0928`、kt1 の `/mnt/data1/ecalj_test0930b`・`0930c`）も trash に入れるか
- ブランチ `fix-idu10`（main にマージ済み）を消すか
- push: dev・rel とも、t14 の main より 859 コミット遅れ（2026-10-02 14:35。2026-10-02 に未公開の範囲の履歴を書き換えた）

### 実行中

- Löwdin を標準にした版（`2ffb0331f`）の全部の試験の組: kr7 `~/ecalj_testL`、kt1 `/mnt/data1/ecalj_testL`（2026-10-02 14:35 の時点で bench の組の途中）。
  済んだ組は両方で PASS（inputs、install、eps、procar、afsym、affix、samples の 13）。mlo・mloqsgw は古い参照と比べて FAIL（予想どおり）→ 終わったら新しい参照で
  mlo・mloqsgw・magnon を回し直す。t14 は新しい参照で mlo 45・mloqsgw 5・magnon 2・install 66 が PASS（`~/work/tests_lowdin3`）

---

## 2. TODO（済）（新しい順）

### 2026-10-02 の TODO の進め方と結果（05:24 の計画、2026-10-02 14:35 に済みとした）

計画（user「TODO の手順はよく考えて順序立てて」「まかせる」）:

依存の向き: 対称性の整理（S1〜S6）が、MLO の最大局在化（Sakuma 型の拘束に操作が要る）と AF の項目の土台。最大局在化の結果で、マグノンの窓と
Wannier とのずれの項目の判断が変わる。だから対称性 → 最大局在化 → マグノン・ずれの順。長い試験の間に、独立した小さな項目を挟む。

| 順 | 項目 | なぜこの位置か | 確かめ | 状況（2026-10-02 12:19） |
| --- | --- | --- | --- | --- |
| 1 | 対称性 S1（分ける、等価変換） | 以後の対称性の段の土台 | 172 入力の lmchk が S0 と同じ、試験の組 | **済み** `1472f3ad7` |
| 2 | 小さい独立の直し: `job_mlo_soc` の空の spin2、`m_tetrakbt` の使われないルーチン、`auto_creplot.py` | 1・3 の試験の待ち時間に | mlo、kBT の組 | 前二つ**済み**。`auto_creplot.py` は移すか捨てるか要判断 |
| 3 | 対称性 S2（GW 側の `mptauof` の重複を外す） | S3 の前に使う側を一本に | gwall がビット単位で同じ | **済み**（最小の形）`7f725713b` |
| 4 | GPU の build の module の循環 | kt1・kr7 で回せる、手元と並行 | kt1・kr7 の clean build と試験 | **済み** `dad040932`・`8d8f7b880` |
| 5 | 対称性 S3（`symmetry.json` を読む口）→ S4（純粋な並進）→ S5（AF、`AFsymmetry/NiO` の pwmode も）→ S6（既定に） | 設計どおり一段ずつ | 各段の表（`MD/symmetry_spglib.md` §4.7） | S3〜S6 **済み**（S6: lmf が同梱の spglib で求める、user の判断 2026-10-02。2026-10-02 08:34 から 3 台で試験） |
| 6 | MLO の最大局在化（Python で試作、Fe・Ni） | 5 の操作を使う | Ω、U、マグノン | 試作**済み**（`mlo_maxloc.py --sym`）。Ω まで。U・マグノンは未 |
| 7 | マグノンの既定の窓、MLO と Wannier のずれ | 6 の結果で判断 | Fe・Ni・FeCo | **済み**: 窓は既定 (2, 2) のまま。Löwdin にすると窓によらず Wannier 版に合う（研究ログ 11:12）。Löwdin を MLO の標準にした（`d10a63716`、研究ログ 12:16） |
| 8 | MLO の模型の残り（EH2 の崩れ、§9 の目安の値の測り直し、空格子球の自動化） | 計算機で裏で回せる | MATERIALS | 原因の確かめと測り直し、`ctrlg_addes.py` の直し**済み**。EH2 の直し方は要判断 |
| 後 | `sugw` のメモリ、`hgw` の残り、MP の API と GW1500 の選定 | 大きい、または外の事情 | — | 未 |

- 済んだ項目（2026-10-02 14:35 に §1 から移した）: Löwdin の後の文書（ecaljdoc mlo §6・表 M8、Changes.txt (6)、handover、`Fe_mlo_magnon/README.md`）。AFTEST の (d)（使い方を ecaljdoc の UsageDetailed.md へ）。
  マグノンの既定の窓（窓は (2, 2) のまま。Löwdin で窓によらない）。対称性 S0〜S6（`MD/symmetry_spglib.md` §4.7）。kt1 の GW1500 の run3（2026-10-01 16:43 に終了）

### 2026-10-02（やったこと）

- user の判断（2026-10-02）を実行した（2026-10-02 12:59）: (1) 旧 `a0c7a7300`（新 `6e2903731`）に誤って入っていたビルドの生成物 3197 本を、未公開の範囲の書き換えで履歴から除いた
  （main のツリーは同じ、dev・rel のコミットは同じ。文書とコミットメッセージのハッシュは直した、対応表 `MD/commit_map_20261002.txt`、控えは TAKAOMINI）。`.git` 1.1 GB → 776 MB。
  (2) ecaljdoc の古い文書（`BackUp/`、`ecaljdetails/`、古い書き出しの pdf）を ecaljdoc の trash へ、要点は past_log.md §14。(3) trash は適宜減らす（`ecaljclaude.md`）。
  API キーは今のツリーと `cb9b2d7b7` 以後のコミットに無いことを確かめた
- **Löwdin で直交化した MLO を標準に**（2026-10-02 12:19、user「Löwdin 直交化を MLO の標準に（バンドは変わらない）」「全体的に調べて、O なしで OK ならそっちをメイン」
  「やれた範囲での決断として O なし」）: 模型は H̃(R) だけ（O(R) = δ）、__cmlo・sugw の a'・m_sigmlo も同じ基底（`d10a63716`）。63 物質で E_F 近くの最大のずれの中央値
  0.019 → 0.011 eV、PASS 53 → 55（研究ログ 12:16、表 12:16-1）。読み込み時の `--mlo_lowdin`（`d46c4014d`、11:12 のマグノン・cRPA）は消した。前の非直交の模型は `--mlo_raw`。
  §2 の質問「マグノン・U の既定を Löwdin にするか」はこれで閉じた
- 対称性 S6（2026-10-02 08:34、user「lmf でつくればいい」）: spglib 2.6.0 の C を同梱し、lmf・lmchk が操作を求めて `symmetry.<sname>.json` を書く（`67f9c1b03`）。172 入力で Python 版と同じ
- ecaljdoc mlo §9 の式 (12) の目安の値を測り直した（06:55、`~/work/ovlp_20261002`、63 物質の MLO の段だけ）: 規格化の後は最小 0.20〜0.42、中央値 0.24〜0.46。FAIL の 10 物質は前と同じ（ecaljdoc `mlo.md`、TODO から外した）
- 対称性（`MD/symmetry_spglib.md` §4.7）: S0 spglib との照合（`b889a00b8`、172 入力で食い違い 0）、S1 分割（`1472f3ad7`）、S2（`7f725713b`）、
  S3 `symmetry.<sname>.json` を読む口（`7902490e4`）、S4a `mptauof` が渡された並進を使う（`6ff09964e`）。S4b（純粋な並進を操作に）は試験中
- GPU の build の module の循環を切った（`dad040932`）、CMake の回避策を外した（`8d8f7b880`）。kt1・kr7 でまっさらなビルドが通った（TODO から外した）
- `m_tetrakbt` の使われないルーチン（`4e2ea8957`）、SOC の MLO が空の spin2 を書く件（`65a5c903f`）（TODO から外した）
- MLO を実空間で規格化（`cbb81dede`）、`mlo_spread.py` と比較の記録（`6e0520052`）、比較の一式とタグ `last-wannier`（`dbcd6e51d`）、cRPA の試験を MLO 版に（`17ac5f104`）、Wannier・AHC・lmfham2 を外した（`5e7244eff`、`4d0989151`）。`MD/wannier_vs_mlo.md`、ecaljdoc mlo §6
- MLOsamples などの古い作業ファイル 68 本（`ctrlp.*`、`lmfham2parameters.check`、`out_lmfham1`、`bandplot_MPO.*` など）を trash へ（past_log 表 1）
- `TOOLS/sync_ecalj_src.sh`: 送り先の `SRC/subroutines`・`main`・`exec` にあって HEAD に無いファイルを `trash/` へ移す（Wannier を外したとき kt1・kr7 に 30 本残って CMake が拾った）

### 2026-10-01

- Samples の旧形式の入力を trash へ（user の指示）: TestInstall の `ctrl.*`（`d436935ce`）、残りの Samples の `ctrl.*` 28 本（`b54a11003`）、MATERIALS の `GWinput`・`ctrlgenM1.ctrl.batio3` と `MLOsamples/Al2O3_Cr/CASE1ok`〜`CASE5ok`（`5c954efcf`）。中身の要点は past_log.md 表 1
- MLO の模型の検査（`mlo_bandcheck.py` の CHECK、`job_mlo`・`job_mlo_soc` が最後に回す、`mlo` が `MLO_ovlpmin.dat` を書く）、
  窓の基準の E_F を `--efermi=` のファイルから取る直し、`Samples/MATERIALS` の 63 の ctrlg の `[mlo]` を今の gwinit で書き直し（研究ログ 2026-10-01 21 時）。
  試験: mlo 45、MLO-QSGW 5、install 64、inputs 176 が t14 で PASSED
- MLO の模型を基準 1・2・3 に整理した（user と決めた。ecaljdoc mlo §1・§4・§9、`MD/handover.md` §5、研究ログ 2026-10-01 16〜20 時）:
  半内殻の局所軌道は帯の上端で自動（E_F − 8 eV より上は EH と入れ替え、−17〜−8 eV は加える、`1ef5c7a68`）、gwinit は f と `!` 付きの `mlo_lm2`（基準 2）を書く、
  MLO-QSGW の凍結した模型の本数の誤り（`5daa42f42`）、`mlo_bandcheck.py` の窓を [VBM − 8, CBM + mlo_delta] に（`760d22467`）、
  `job_mlo_soc` がスピン軌道ありの DFT のバンドも描く（GaAsSoc の −0.1 eV は比べ方の誤り、`ab9570ead`）。Samples/MATERIALS の結果のページ
  https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH（`Samples/MATERIALS/mlocheck/gallery.py`）。試験: inputs 176、mlo 25 試料、MLO-QSGW、Fe_mlo_magnon が t14 で PASSED
- AFTEST の修正（`aftest-fix`）を main にマージ（`aa4c24129`、user「マージできるよね」→「良い方を壊さないように」）。AF の対称性ありの計算は ehf・sev・モーメント・場がマージ前と同じで、ehk だけが制約の下の全エネルギーになった。afsym 4・install 64 が PASS。例を `Samples/AFfixMMOM`（user の命名、`NiO_afsym`・`NiO_noafsym`、試験の組 `affix` 12 件 PASS）に置いた（`aaa944368`）
- ecaljdoc の TOML の流れの最短の手順は、`manual/README_tutorial.md` の GetStarted（Step 0〜6: POSCAR → ctrls → `ctrlgenToml.py` → `lmfa`・`lmf` → バンド → `gwsc`、Step 2-Migration、`--ctrlg:` の上書き）に既にあった。TODO から外した。ecaljdoc `manual/mlo.md` に既定の模型が外れる 2 つの型と表 M1 を書いた（`422dc92`、未 push）
- `Samples/MATERIALS` の 65 物質で LDA と MLO の自動の模型を回した（InAs/GaSb n10 の 40 原子は user の判断で回さない）。結果の表は `Samples/MATERIALS/MLOcheck_20261001*.tsv`、一覧のページ https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH、研究ログ 2026-10-01 朝 07:16
- `job_mlo`: ctrlg の `[ham] so = 1` なら `job_mlo_soc` を案内して止まる（`mlo` がハミルトニアンの NaN で止まっていた。GaAs_so）。`--ctrlg:ham.so=` の上書きは尊重。`mlo_bandplot.py` は空の spin2 を描かない
- `--cls`: `m_clsmode_finalize` に `ndimh`（`m_igv2x` の最後の k 点の値）でなく `nbandmx` を渡す（`vcdmel` の重みの並びと `dostet` の読み方を揃える）。CrN（`Samples/TestInstall/crn`）で APW なし・pwemax 2 と 5・4×4×4（ndimh が 182〜188 と変わる）のどれも `dos-vcdmel.crn` が一致（ずれるのは DOS の窓より上の帯だけだった）。ブランチ `cls-nbandmx` をマージ（`f2e298ac3`）。試験は `~/work/clstest_20261001`
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
  09-30 の idu の修正（`f14270e83`）から lmf で止まっていた（`LDA+U must be spin-polarized!`）ので同じく直した（`66d3edbee`）
- `Samples/MATERIALS` を仕分けた: 構造のデータベース（62 物質）を `Database/` に展開（`ctrls` と、GW・MLO の節つきの `ctrlg`、`lmchk` で全部読める）。ほかの旧形式の 30 項目は trash（624 ファイル）。残したのは `La2CuO4`・`InAsGaSb`・`BaTiO3`。拾ったノウハウは past_log.md §9
- `README.txt` と `.org` を Markdown に（user「README.txt とあるのは md 形式に。org もそう」）: `StructureTool/README.txt`、`Samples/EPS/EPS_GaAs/README_eps.org`、`Samples/MATERIALS/{LaGaO3_relax,MLOsamples,NiSe_aftest}/README.org`、`Samples/MATERIALS/Si_doping_sample/Memo_bgcharge.org`（`git mv`）。pandoc は `_` を下付きに、`--` をダッシュに変えるので使わず、原文を保つ変換（見出し・`#+TITLE`・`#+begin_src`・コマンド行と設定の断片をコードブロックに）
- `GetSyml/`・`StructureTool/` の古い例（旧形式の `ctrl.*` 86 本、`syml.*`、鉱物の POSCAR 135 本）と使わないスクリプトを trash へ（257 ファイル）。本体（`getsyml`、`vasp2ctrl`・`ctrl2vasp`、`viewvesta`、`refineposcar.py`、`superlattice/`）は残し、動くことを確かめた。past_log.md §12
- `Doxygen/` を trash へ（user「そうしよう」）。コメントを Doxygen 形式に揃えることはしない。作り直し方は past_log.md §11
- `TOOLS/` の古い道具（約 60 項目、632 ファイル）を trash へ。残したのは `samples_tests.sh`・`sync_ecalj_src.sh`・`ozbench/`。中身は past_log.md §10。`diffnum` は試験で今も使うが、使うのは `SRC/exec/pylib/diffnum0.py`（TOOLS の版は古い）
- Claude の個人メモリから、引き継ぐ価値があり今も正しいものを [handover.md](handover.md) に写した（user「メモリの内容はパッケージに入らないので MD/ に。混乱を招くものは良くない」）。
  写す前に 5 点をコードと照らした（rel の既定ブランチは `main`、`master` は無い、など）
- `README.md` を `MD/README.md` へ（user「Claude に読ませて、人間は Claude から情報を取る構造にする。README を人間に読ませるのは好ましくない」）。
  最上位の `README.md` は数行の案内だけ（GitHub の表紙が空にならないように）。`ecaljclaude.md` の方針の文言を直した
- 開発の文書を `MD/` に（`ecaljclaude.md`、`TODOandQuestion.md`、`past_log.md`、研究ログ `research_log.md`、`d72216e17`）
- ecalj の片付けの 3 回目: `PHASE1B_REFACTOR.md`、`HIGHLIGHTS_2026-06_09.md`、`FiniteT_and_QPE_HOWTO.md`、`ecaljdoc_drafts/`、`jobauto/`、
  `SRC/exec_legacy/` を trash へ。中身は past_log.md §3.2・§4.3・§5〜§8 に整理し、README と ecaljdoc（ForDevelopers・kBT・README_tutorial）の参照を直した。
  最上位の `MATERIALS/` は trash ではなく `Samples/MATERIALS/` の下へ移した（user「いったん Samples の下へ」、624 ファイル、`git mv`）
- ecalj の片付け（user「いらないものは trash へ。ノウハウは過去ログへ。重複は整理してから trash へ」）:
  ビルドの生成物 1853 ファイル（`df4968a5e`）、最上位の打ち込み用と古いスクリプト、`.refactor_notes/`、`SRC/exec/BK`、
  MLOsamples の古い試行、`TestInstall/TESTunused`、`*.bk` などを trash へ。ノウハウは [past_log.md](past_log.md)
- 最上位の `SRC/BK`・`SRC/execgfortran`・`SRC/execAHC`（181 MB）と、GW1500 の古いスクリプトを追跡から外して trash へ（`0c8dbaaf2`）
- kt1 の GW1500 の作業ファイルと退避（255 GB）を `~/ecalj/trash` へ（`47ce04cc5`）
- ecalj・ecaljdoc の控えのクローンを外付けドライブ `/media/takao/TAKAOMINI`（`ecalj_20261001`・`ecaljdoc_20261001`）に
- Materials Project の API キーを追跡から外した。`<ecalj>/MaterialProject.key`（無視される）と見本 `.example`、`pylib/mpkey.py`（`cb9b2d7b7`）
- GW1500: 物質ごとの注記の表 `ecalj_auto/gw1500_notes_20261001.tsv`、分類と条件の違い（`GW1500_status.md` §5、`b0fdfda19`）。
  5 月に収束した GOOD はやり直さない（run3 のキューを NOTCONV だけに）。Rb8（mp-1179832）は不正な構造と分かり止めた（`221e0e307`）
- `gwscconv`: ギャップの無い反復（金属）は固有値の変化で収束を判定（`1730f7e7e`）
- kt1・kr7（nvfortran GPU）・mic（ifx）で試験の組 inputs・mlo・afsym・install・samples がすべて PASS（HEAD `b695fa65c`）

### 2026-09-30

- `fix-idu10` をマージ（`f14270e83`）: `idu = 10 + mode` が `sigm` の無いとき二重計数なしの LDA+U になっていた（2023-09 から）。
  MLOsamples の SmP・GdCo5 の LDA+U を作り直した（`6e5817ffe`、`23429f0a3`）
- GW1500: FAILED の 147 物質を fp32 で回し直し、144 が収束（`70f79b23b`）
- Samples/Legacy を片付けて組み直し（AtomDimer/N2 など）、ecalj と ecaljdoc を clone し直した
