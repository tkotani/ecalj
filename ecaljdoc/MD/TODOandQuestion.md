# TODO と質問、やったこと

直すべき点は見つけてもその場では直さず、ここに書く（user 2026-10-01）。メンテナに決めてほしいことも、やったことも、ここに書く。
入口は [CLAUDE.md](../../CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）。片付けたものの中のノウハウは [past_log.md](past_log.md)。
細かい経緯と数値は研究ログ [research_log.md](research_log.md)。

書き方: §1 はまだのものだけ（各項目に見つけた日）。済んだら §2 へ移し、済んだ日とコミットを書く（user 2026-10-02「まだのものを冒頭に、終わったものは 2. TODO（済）に」）。
計算機ごとの片付け（trash の削除、ディスクの空き）はここに書かない（user 2026-10-02「ローカル情報はコミットに入れなくていい」）。

---

## 1. TODO（未）（2026-10-02 20:01 に書き直した）

各項目の頭の印: **【やりかけ】** 手を付けて途中、**【未着手】**、**【判断待ち】** user に決めてほしい、**【残す】** user が残すと決めたもの（急がない）。

### MLO

- **【未着手】EuO に EH2 を足すと模型を作る所で止まる**（2026-10-01）: `Hreduction` の規格化の確かめ（band 26、シードのノルムの減り −1.3 %）。EH と EH2 の両方をシードにすると
  `zhev_tk4` が一次従属に近い向きを落とす。Löwdin の模型でも同じ。4f の原子は基準 2 から除く決まり（§9）なので急がない
- **【残す】η ≠ 1 の原因**（2026-10-02）: Goldstone の条件の倍率が Löwdin で Fe 1.24、FeCo 1.28、Ni 1.75。Ni が大きい理由は未確認
- **【残す】空格子球の自動化（基準 3）**（2026-10-01）: [`SRC/exec/ctrlg_addes.py`](../../SRC/exec/ctrlg_addes.py)（`8c158ee96`、像の数を面間隔から決める直しは済み）を ctrlg の生成に組み込むか。
  Bi₂Te₃ の空隙は 2.61 a.u. で 3.0 未満、SiO₂ は手で置いて 0.001 eV
- **【未着手】<u>Fe のマグノン: 小さい q で 2 割高い件</u>と、実験との比較**（2026-10-02）: Löwdin の MLO は q ≤ 0.3 で Wannier 版より 2 割高い（q = 0.1 で 0.085 対 0.068 eV）。
  q ≥ 0.4 は合う。見る所: [`Samples/Magnon/Fe_mlo_magnon/README.md`](../../Samples/Magnon/Fe_mlo_magnon/README.md)（表 2・図 1、`magnon_peaks.npz`）、[`ecaljdoc/MD/wannier_vs_mlo.md`](wannier_vs_mlo.md) §3a（表 3a）、研究ログ 2026-10-02 11:12・
  09:42・07:10・04:17（窓と Löwdin）、作業の場所 `~/work/magnon_nifeco/magnon_window.md`（Fe・FeCo・Ni の窓と Löwdin の図、ページ https://claude.ai/artifact/A9r4CEPHgD1UkYsVA5XQKH）。
  実験の分散（中性子散乱）のデータはリポジトリに無い。文献から取る

### 対称性

- **【やりかけ】操作の順番で QSGW の結果が変わる**（2026-10-02）: heavy の `nio_gwsc444`（AF NiO、R-3m、操作 12）で、spglib（今の既定）と従来の探し方は同じ 12 個の操作を
  違う順番で並べる（従来 e, i·r3d, r3, i, …、spglib e, i, r3d, i·r3d, …）。既約な k・重み・四面体の数・LDA は一致するが、QSGW の 1 反復の QP が最大 15 meV 違う
  （kr7・kt1 で同じ。従来の探し方なら参照と 3 meV）。有限温度の四面体法ではない（この試験は χ0 が T = 0）。物理は順番によらないはずなので、GW のどこかが
  操作の選び方に依っている（候補: q を代表へ移す回転に最初の操作を使う所、縮退や `emax_sigm` の切れ目、オフセット Γ の点）。どこかは未特定。
  **【判断待ち】** 特定するまで、参照を新しい値にするか

### 有限温度

- **【未着手】有限温度の四面体法（`t_tetrakbt > 0`、`m_tetrakbt`）と従来の T = 0 の四面体法（`t_tetrakbt = 0`、`tetwt5`）の関係を確かめる**（2026-10-02、user）:
  T → 0 で有限温度の重みが従来の重みに一致するか（同じ k メッシュ・同じ対称性で、χ0 の虚部と QP のエネルギーを比べる。T を 300、100、30、10 K と下げる）。
  Σ 側の Fermi–Dirac の幅（`t_sigmaw`）と χ0 側の温度の組み合わせも

### コード

- **【未着手】HEAD に残っている生成物らしいもの**（2026-10-02、旧 `a0c7a7300` で入り今も追跡）: `SRC/.#Memo4rotation`、[`SRC/exec/cmake_install.cmake`](../../SRC/exec/cmake_install.cmake)、`hello.py`、`platform`（800 KB）、
  `lmf2.py`・`lmchk.py`・`pylmfa`・`pysample`・`ohtaka`・`epsPPd`・`epsPPsaito`・`job_senefbz`・`readeps_dig2.py`・`auto_kauto.py`。使われているかを確かめて trash へ
- **【未着手】[`SRC/exec/auto_creplot.py`](../../SRC/exec/auto_creplot.py)**（2026-10-01）: 旧形式の `ctrl.<sname>` を書き換える。`auto_job_mp.py` が使う。ctrlg に直すか、使わないなら trash へ
- **【未着手】`sugw`（`lmf --jobgw=1`）のメモリ**（2026-10-01）: `GEIGpart` の `ppovl(ngp,ngp)` と `ppovlLU` を各ランクで持つ（32·ngp² バイト、胞の体積の 2 乗）。
  案: (a) ngp から並列数を決める、(b) Cholesky の因子だけ持つ、(c) O·x を FFT で作り反復法で解く
- **【未着手】`hgw` の残り**（2026-04、past_log.md §3.2）: ノード内の W の共有（`MPI_Win_allocate_shared`）、`hsfp0_sc` の Sx・core の交換も `hgw` に（優先度は低い）
- **【未着手】AFTEST の残り**（2026-10-01、研究ログ 2026-10-01 06:46）: (a) モーメントを下げる向きで更新が行き過ぎる（割線で見積もるか）、(b) 対がサイト 1・2、ブロック 1・2 の決め打ち、
  (c) `m_ldau_init` が lmf の起動のたびに場を更新して `mmagfield.aftest` を書き直す（`job_band` でも）

### GW1500 と Materials Project

- **【未着手】`auto_mpquery.py` が今の MP で動かない**（2026-10-01）: 新しい ID の形で `mp_api` の検証が止まる。`mp_api` を上げるか、REST を直接読む
- **【未着手】GW1500 の選定に構造の確かめを入れる**（2026-10-01）: 副格子を抜き出した MP の項目が選ばれていた（[`ecalj_auto/GW1500_status.md`](../../ecalj_auto/GW1500_status.md) §5.3 の表 8・9）。`auto_mpquery.py` で弾くか印を付ける

### 試験と入力

- **【未着手】`MLOsamples/RuO2` の `rst`・`dmats`** は `pwmode = 11` の LDA+U の誤りの時期（2026-03-30〜09-30）に作ったもの。作り直すか
- **【判断待ち】試験の入力の温度**（`t_tetrakbt = 262`、`t_sigmaw = 0`）をテンプレート（300/300）に揃えるか。揃えると gas_gwsc・fe_gwsc などの参照が動く

### 判断待ち（ほかに）

- **GW1500 の `INVALID_STRUCTURE`（12）・`SUSPECT_STRUCTURE`（4）を集合から外すか**。いまは注記だけ
- 試験に使った古いツリー（mic の `~/ecalj_test0928`、kt1 の `/mnt/data1/ecalj_test0930b`・`0930c`）も trash に入れるか
- ブランチ `fix-idu10`（main にマージ済み）を消すか
- push: dev・rel とも、t14 の main より 888 コミット遅れ（2026-10-02 21:07。2026-10-02 に未公開の範囲の履歴を書き換えた）。文書の公開（ecalj/ecaljdoc への subtree push）も同じとき

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
| 5 | 対称性 S3（`symmetry.json` を読む口）→ S4（純粋な並進）→ S5（AF、`AFsymmetry/NiO` の pwmode も）→ S6（既定に） | 設計どおり一段ずつ | 各段の表（[`ecaljdoc/MD/symmetry_spglib.md`](symmetry_spglib.md) §4.7） | S3〜S6 **済み**（S6: lmf が同梱の spglib で求める、user の判断 2026-10-02。2026-10-02 08:34 から 3 台で試験） |
| 6 | MLO の最大局在化（Python で試作、Fe・Ni） | 5 の操作を使う | Ω、U、マグノン | 試作**済み**（`mlo_maxloc.py --sym`）。Ω まで。U・マグノンは未 |
| 7 | マグノンの既定の窓、MLO と Wannier のずれ | 6 の結果で判断 | Fe・Ni・FeCo | **済み**: 窓は既定 (2, 2) のまま。Löwdin にすると窓によらず Wannier 版に合う（研究ログ 11:12）。Löwdin を MLO の標準にした（`d10a63716`、研究ログ 12:16） |
| 8 | MLO の模型の残り（EH2 の崩れ、§9 の目安の値の測り直し、空格子球の自動化） | 計算機で裏で回せる | MATERIALS | 原因の確かめと測り直し、`ctrlg_addes.py` の直し**済み**。EH2 の直し方は要判断 |
| 後 | `sugw` のメモリ、`hgw` の残り、MP の API と GW1500 の選定 | 大きい、または外の事情 | — | 未 |

- 最大局在化（MV）の試作は [`TOOLS/gadget/`](../../TOOLS/gadget/README.md)（README つき、bindir に入らない）に移して残した（2026-10-02 16:05、user「メインにしない。役に立ちうるなら TOOLS の下に gadget」）。
  役に立ちうる所: 結合の上の局在関数、Wannier90 との比較、Löwdin の局在の目安
- Löwdin を標準にした後の内挿のずれは許容範囲とした（2026-10-02 16:13、user「許容範囲というべき。Löwdin でいい」）。窓 [VBM − 3, CBM + 2] eV で 63 物質の最大のずれの中央値
  0.035 → 0.027、rms 0.0067 → 0.0051、0.1 eV を超えるもの 10 → 8（どれも生の模型でも超える）。悪くなったのは 2H-SiC 0.055 → 0.072、GaAs 0.044 → 0.056、
  GaSb 0.078 → 0.091、MnO 0.010 → 0.030、NiO 0.011 → 0.025。Fe の 4s 帯の底（−8〜−3 eV）は窓の外（`~/work/lowdin_20261002/win3_2.tsv`）
- [`Samples/AFsymmetry/NiO`](../../Samples/AFsymmetry/NiO) を `pwmode = 11` に、AF の入力を `symgrpaf = "find"`（spglib が AF の操作を見つける）に（2026-10-02 16:39、`6f3642d8f`、user「symgrpAF も spglib に見つけさせたい」）。
  従来の探し方で `"find"` は理由を言って止まる。afsym 4・affix 12・mlo 45 が PASS
- `mlo_bandcheck.py` の判定 (3) を「帯のとげ」（模型のバンドの 2 階差分、窓 [VBM − 3, CBM + 2] で 2 eV 超）に替えた（2026-10-02 20:01、user「線形独立が壊れるとき固有値が急に変わる」）。
  壊れた Cu・Ni（EH2、生の模型）は 24・926 eV、健全な 63 物質は最大 0.56 eV。ecaljdoc mlo §9 式 (12)
- 基準 2・3 を Löwdin の模型で確かめた（2026-10-02、`~/work/lowdin_crit23_20261002`）: 生の模型で壊れた Cu・Ni の EH2 は PASS（0.045、0.055 eV）、基準 1 で外れた 7 物質は基準 2 で PASS、
  SiO₂ は基準 3 で 0.005 eV。EH2 の崩れの直し方（正準直交化など）は要らなくなった（EuO だけ §1 に残る）
- MLO と Wannier の cRPA の U の差（Ni の d、Löwdin 2.90 対 Wannier 3.78 eV）は、部分空間の取り方の違いで、直すものではないとした（2026-10-02 20:13、user）。RPA の U は
  基底にほとんど依らない（1.43 対 1.58）。[`ecaljdoc/MD/wannier_vs_mlo.md`](wannier_vs_mlo.md) §3a、ecaljdoc mlo §6
- Löwdin を標準にした版の試験（2026-10-02 21:07）: kr7・kt1 で全部の組が PASS（inputs 172、install 66、eps 18、procar 5、afsym 4、affix 12、samples の 13、bench 2、
  新しい参照で mlo 45・mloqsgw 5・magnon 2）。heavy の `nio_gwsc444` だけ 15 meV（§1「対称性」）
- 済んだ項目（2026-10-02 14:35 に §1 から移した）: Löwdin の後の文書（ecaljdoc mlo §6・表 M8、Changes.txt (6)、handover、`Fe_mlo_magnon/README.md`）。AFTEST の (d)（使い方を ecaljdoc の UsageDetailed.md へ）。
  マグノンの既定の窓（窓は (2, 2) のまま。Löwdin で窓によらない）。対称性 S0〜S6（[`ecaljdoc/MD/symmetry_spglib.md`](symmetry_spglib.md) §4.7）。kt1 の GW1500 の run3（2026-10-01 16:43 に終了）

### 2026-10-02（やったこと）

- user の判断（2026-10-02）を実行した（2026-10-02 12:59）: (1) 旧 `a0c7a7300`（新 `6e2903731`）に誤って入っていたビルドの生成物 3197 本を、未公開の範囲の書き換えで履歴から除いた
  （main のツリーは同じ、dev・rel のコミットは同じ。文書とコミットメッセージのハッシュは直した、対応表 [`ecaljdoc/MD/commit_map_20261002.txt`](commit_map_20261002.txt)、控えは TAKAOMINI）。`.git` 1.1 GB → 776 MB。
  (2) ecaljdoc の古い文書（`BackUp/`、`ecaljdetails/`、古い書き出しの pdf）を ecaljdoc の trash へ、要点は past_log.md §14。(3) trash は適宜減らす（[`ecaljclaude.md`](ecaljclaude.md)）。
  API キーは今のツリーと `cb9b2d7b7` 以後のコミットに無いことを確かめた
- **Löwdin で直交化した MLO を標準に**（2026-10-02 12:19、user「Löwdin 直交化を MLO の標準に（バンドは変わらない）」「全体的に調べて、O なしで OK ならそっちをメイン」
  「やれた範囲での決断として O なし」）: 模型は H̃(R) だけ（O(R) = δ）、__cmlo・sugw の a'・m_sigmlo も同じ基底（`d10a63716`）。63 物質で E_F 近くの最大のずれの中央値
  0.019 → 0.011 eV、PASS 53 → 55（研究ログ 12:16、表 12:16-1）。読み込み時の `--mlo_lowdin`（`d46c4014d`、11:12 のマグノン・cRPA）は消した。前の非直交の模型は `--mlo_raw`。
  §2 の質問「マグノン・U の既定を Löwdin にするか」はこれで閉じた
- 対称性 S6（2026-10-02 08:34、user「lmf でつくればいい」）: spglib 2.6.0 の C を同梱し、lmf・lmchk が操作を求めて `symmetry.<sname>.json` を書く（`67f9c1b03`）。172 入力で Python 版と同じ
- ecaljdoc mlo §9 の式 (12) の目安の値を測り直した（06:55、`~/work/ovlp_20261002`、63 物質の MLO の段だけ）: 規格化の後は最小 0.20〜0.42、中央値 0.24〜0.46。FAIL の 10 物質は前と同じ（ecaljdoc `mlo.md`、TODO から外した）
- 対称性（[`ecaljdoc/MD/symmetry_spglib.md`](symmetry_spglib.md) §4.7）: S0 spglib との照合（`b889a00b8`、172 入力で食い違い 0）、S1 分割（`1472f3ad7`）、S2（`7f725713b`）、
  S3 `symmetry.<sname>.json` を読む口（`7902490e4`）、S4a `mptauof` が渡された並進を使う（`6ff09964e`）。S4b（純粋な並進を操作に）は試験中
- GPU の build の module の循環を切った（`dad040932`）、CMake の回避策を外した（`8d8f7b880`）。kt1・kr7 でまっさらなビルドが通った（TODO から外した）
- `m_tetrakbt` の使われないルーチン（`4e2ea8957`）、SOC の MLO が空の spin2 を書く件（`65a5c903f`）（TODO から外した）
- MLO を実空間で規格化（`cbb81dede`）、`mlo_spread.py` と比較の記録（`6e0520052`）、比較の一式とタグ `last-wannier`（`dbcd6e51d`）、cRPA の試験を MLO 版に（`17ac5f104`）、Wannier・AHC・lmfham2 を外した（`5e7244eff`、`4d0989151`）。[`ecaljdoc/MD/wannier_vs_mlo.md`](wannier_vs_mlo.md)、ecaljdoc mlo §6
- MLOsamples などの古い作業ファイル 68 本（`ctrlp.*`、`lmfham2parameters.check`、`out_lmfham1`、`bandplot_MPO.*` など）を trash へ（past_log 表 1）
- [`TOOLS/sync_ecalj_src.sh`](../../TOOLS/sync_ecalj_src.sh): 送り先の [`SRC/subroutines`](../../SRC/subroutines)・`main`・`exec` にあって HEAD に無いファイルを `trash/` へ移す（Wannier を外したとき kt1・kr7 に 30 本残って CMake が拾った）

### 2026-10-01

- Samples の旧形式の入力を trash へ（user の指示）: TestInstall の `ctrl.*`（`d436935ce`）、残りの Samples の `ctrl.*` 28 本（`b54a11003`）、MATERIALS の `GWinput`・`ctrlgenM1.ctrl.batio3` と `MLOsamples/Al2O3_Cr/CASE1ok`〜`CASE5ok`（`5c954efcf`）。中身の要点は past_log.md 表 1
- MLO の模型の検査（`mlo_bandcheck.py` の CHECK、`job_mlo`・`job_mlo_soc` が最後に回す、`mlo` が `MLO_ovlpmin.dat` を書く）、
  窓の基準の E_F を `--efermi=` のファイルから取る直し、[`Samples/MATERIALS`](../../Samples/MATERIALS/README.md) の 63 の ctrlg の `[mlo]` を今の gwinit で書き直し（研究ログ 2026-10-01 21 時）。
  試験: mlo 45、MLO-QSGW 5、install 64、inputs 176 が t14 で PASSED
- MLO の模型を基準 1・2・3 に整理した（user と決めた。ecaljdoc mlo §1・§4・§9、[`ecaljdoc/MD/handover.md`](handover.md) §5、研究ログ 2026-10-01 16〜20 時）:
  半内殻の局所軌道は帯の上端で自動（E_F − 8 eV より上は EH と入れ替え、−17〜−8 eV は加える、`1ef5c7a68`）、gwinit は f と `!` 付きの `mlo_lm2`（基準 2）を書く、
  MLO-QSGW の凍結した模型の本数の誤り（`5daa42f42`）、`mlo_bandcheck.py` の窓を [VBM − 8, CBM + mlo_delta] に（`760d22467`）、
  `job_mlo_soc` がスピン軌道ありの DFT のバンドも描く（GaAsSoc の −0.1 eV は比べ方の誤り、`ab9570ead`）。Samples/MATERIALS の結果のページ
  https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH（[`Samples/MATERIALS/mlocheck/gallery.py`](../../Samples/MATERIALS/mlocheck/gallery.py)）。試験: inputs 176、mlo 25 試料、MLO-QSGW、Fe_mlo_magnon が t14 で PASSED
- AFTEST の修正（`aftest-fix`）を main にマージ（`aa4c24129`、user「マージできるよね」→「良い方を壊さないように」）。AF の対称性ありの計算は ehf・sev・モーメント・場がマージ前と同じで、ehk だけが制約の下の全エネルギーになった。afsym 4・install 64 が PASS。例を [`Samples/AFfixMMOM`](../../Samples/AFfixMMOM/README.md)（user の命名、`NiO_afsym`・`NiO_noafsym`、試験の組 `affix` 12 件 PASS）に置いた（`aaa944368`）
- ecaljdoc の TOML の流れの最短の手順は、`manual/README_tutorial.md` の GetStarted（Step 0〜6: POSCAR → ctrls → `ctrlgenToml.py` → `lmfa`・`lmf` → バンド → `gwsc`、Step 2-Migration、`--ctrlg:` の上書き）に既にあった。TODO から外した。ecaljdoc `manual/mlo.md` に既定の模型が外れる 2 つの型と表 M1 を書いた（`422dc92`、未 push）
- [`Samples/MATERIALS`](../../Samples/MATERIALS/README.md) の 65 物質で LDA と MLO の自動の模型を回した（InAs/GaSb n10 の 40 原子は user の判断で回さない）。結果の表は `Samples/MATERIALS/MLOcheck_20261001*.tsv`、一覧のページ https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH、研究ログ 2026-10-01 朝 07:16
- `job_mlo`: ctrlg の `[ham] so = 1` なら `job_mlo_soc` を案内して止まる（`mlo` がハミルトニアンの NaN で止まっていた。GaAs_so）。`--ctrlg:ham.so=` の上書きは尊重。`mlo_bandplot.py` は空の spin2 を描かない
- `--cls`: `m_clsmode_finalize` に `ndimh`（`m_igv2x` の最後の k 点の値）でなく `nbandmx` を渡す（`vcdmel` の重みの並びと `dostet` の読み方を揃える）。CrN（[`Samples/TestInstall/crn`](../../Samples/TestInstall/crn)）で APW なし・pwemax 2 と 5・4×4×4（ndimh が 182〜188 と変わる）のどれも `dos-vcdmel.crn` が一致（ずれるのは DOS の窓より上の帯だけだった）。ブランチ `cls-nbandmx` をマージ（`f2e298ac3`）。試験は `~/work/clstest_20261001`
- `InstallAll.py`: build の後に、bindir の中で「この ecalj の木を指していて行き先の無いリンク」だけを消す（`remove_dangling_links`）。[`SRC/exec`](../../SRC/exec) の行き先の無いリンク（消した `SRC/exec/build/` を指す 12 本）は張らず、trash に移した。これが `~/bin/mlo` などを一時的に死んだリンクで上書きしていた（CMake の deliver が build の後に戻していた）。SRC/exec のエディタの一時ファイル（`~` で終わる、`#`・`.#` で始まる）はリンクしない。t14 の `~/bin` の 58 本は次のインストールで消える（一時の bindir で試験）
- ecaljdoc `manual/spectrum.md` に、メッシュの外の q（`QforGW`）の注意 3 点（`EMAXforGW` が必須、窓を変えたら交換から、`epsWVR` の行）を書いた（past_log.md §5 から、ecaljdoc の未 push のコミット）
- AFTEST を調べた（研究ログ 2026-10-01 朝 06:46）: afsym ではモーメントを目標に保てるが、表示の ehk に −uhx·m_d、ehf に −2·uhx·m_d が残る。afsym なしでは誤り（サイトの電荷が分かれる）。修正はブランチ `aftest-fix`
- GW1500: kr7 で fp32 と TF32 を同じバイナリ・同じ入力で比べた（8 物質、02:03〜06:24）。最終のギャップの差は 1 meV 未満。5 月との 0.29〜1.48 eV の差は精度ではなく、5 月の振動と設定の違い（`GW1500_status.md` §5.1 の表 7）。GOOD 1120 を精度の理由で見直す必要は無い
- `SRC/subroutines/m_tetrakbt_BUGREPORT.md`（2026-06、直し済みの不具合の報告）の要点を past_log.md §13 に移して trash へ。空の道標 `SRC/TestInstall_is_moved_to_under_ecaljSamples` も trash へ
- `.gitignore` に `/build/`（VSCode の CMake 拡張が最上位に作る）
- [`ecaljdoc/MD/module_map.md`](module_map.md)（生成物）と [`TOOLS/module_map.py`](../../TOOLS/module_map.py): 主プログラム → 入口の module、module の階層（226 module、最大 29 段、`m_lmf` が頂点）、
  依存と被依存の数、module の冒頭のコメント。[`ecaljclaude.md`](ecaljclaude.md) から参照。GPU の build の module の循環が 1 つ見つかった（上の TODO）
- [`GetSyml/README.md`](../../GetSyml/README.md)・[`StructureTool/README.md`](../../StructureTool/README.md)・[`Samples/EPS/EPS_GaAs/README_eps.md`](../../Samples/EPS/EPS_GaAs/README_eps.md) を今の形に書き直した（入力は `ctrlg.<sname>.toml`、
  `getsyml --nobzview`、StructureTool の各スクリプトの向きと出力のファイル名、EPS は `[gw] QforEPS`・`QforEPSau`・`n1n2n3`・`[product_basis] pb_lcutmx`、
  EPS の出力の列）。`getsyml` の使い方の表示も（`-nobzview` → `--nobzview`、`ctrl.nio` → `ctrlg.nio.toml`）
- `ctrlgenToml.py`: nspin=1 のとき原子表の IDU/UH/JH をコメントにして書く（LDA+U は nspin=2 が要る）。`Database/Ce`（nspin=1、idu=12）が
  09-30 の idu の修正（`f14270e83`）から lmf で止まっていた（`LDA+U must be spin-polarized!`）ので同じく直した（`66d3edbee`）
- [`Samples/MATERIALS`](../../Samples/MATERIALS/README.md) を仕分けた: 構造のデータベース（62 物質）を `Database/` に展開（`ctrls` と、GW・MLO の節つきの `ctrlg`、`lmchk` で全部読める）。ほかの旧形式の 30 項目は trash（624 ファイル）。残したのは `La2CuO4`・`InAsGaSb`・`BaTiO3`。拾ったノウハウは past_log.md §9
- `README.txt` と `.org` を Markdown に（user「README.txt とあるのは md 形式に。org もそう」）: `StructureTool/README.txt`、`Samples/EPS/EPS_GaAs/README_eps.org`、`Samples/MATERIALS/{LaGaO3_relax,MLOsamples,NiSe_aftest}/README.org`、`Samples/MATERIALS/Si_doping_sample/Memo_bgcharge.org`（`git mv`）。pandoc は `_` を下付きに、`--` をダッシュに変えるので使わず、原文を保つ変換（見出し・`#+TITLE`・`#+begin_src`・コマンド行と設定の断片をコードブロックに）
- `GetSyml/`・`StructureTool/` の古い例（旧形式の `ctrl.*` 86 本、`syml.*`、鉱物の POSCAR 135 本）と使わないスクリプトを trash へ（257 ファイル）。本体（`getsyml`、`vasp2ctrl`・`ctrl2vasp`、`viewvesta`、`refineposcar.py`、`superlattice/`）は残し、動くことを確かめた。past_log.md §12
- `Doxygen/` を trash へ（user「そうしよう」）。コメントを Doxygen 形式に揃えることはしない。作り直し方は past_log.md §11
- `TOOLS/` の古い道具（約 60 項目、632 ファイル）を trash へ。残したのは `samples_tests.sh`・`sync_ecalj_src.sh`・`ozbench/`。中身は past_log.md §10。`diffnum` は試験で今も使うが、使うのは [`SRC/exec/pylib/diffnum0.py`](../../SRC/exec/pylib/diffnum0.py)（TOOLS の版は古い）
- Claude の個人メモリから、引き継ぐ価値があり今も正しいものを [handover.md](handover.md) に写した（user「メモリの内容はパッケージに入らないので ecaljdoc/MD/ に。混乱を招くものは良くない」）。
  写す前に 5 点をコードと照らした（rel の既定ブランチは `main`、`master` は無い、など）
- [`README.md`](../../README.md) を `ecaljdoc/MD/README.md` へ（user「Claude に読ませて、人間は Claude から情報を取る構造にする。README を人間に読ませるのは好ましくない」）。
  最上位の [`README.md`](../../README.md) は数行の案内だけ（GitHub の表紙が空にならないように）。[`ecaljclaude.md`](ecaljclaude.md) の方針の文言を直した
- 開発の文書を [`ecaljdoc/MD/`](.) に（[`ecaljclaude.md`](ecaljclaude.md)、[`TODOandQuestion.md`](TODOandQuestion.md)、[`past_log.md`](past_log.md)、研究ログ [`research_log.md`](research_log.md)、`d72216e17`）
- ecalj の片付けの 3 回目: `PHASE1B_REFACTOR.md`、`HIGHLIGHTS_2026-06_09.md`、`FiniteT_and_QPE_HOWTO.md`、`ecaljdoc_drafts/`、`jobauto/`、
  `SRC/exec_legacy/` を trash へ。中身は past_log.md §3.2・§4.3・§5〜§8 に整理し、README と ecaljdoc（ForDevelopers・kBT・README_tutorial）の参照を直した。
  最上位の `MATERIALS/` は trash ではなく [`Samples/MATERIALS/`](../../Samples/MATERIALS/README.md) の下へ移した（user「いったん Samples の下へ」、624 ファイル、`git mv`）
- ecalj の片付け（user「いらないものは trash へ。ノウハウは過去ログへ。重複は整理してから trash へ」）:
  ビルドの生成物 1853 ファイル（`df4968a5e`）、最上位の打ち込み用と古いスクリプト、`.refactor_notes/`、`SRC/exec/BK`、
  MLOsamples の古い試行、`TestInstall/TESTunused`、`*.bk` などを trash へ。ノウハウは [past_log.md](past_log.md)
- 最上位の `SRC/BK`・`SRC/execgfortran`・`SRC/execAHC`（181 MB）と、GW1500 の古いスクリプトを追跡から外して trash へ（`0c8dbaaf2`）
- kt1 の GW1500 の作業ファイルと退避（255 GB）を `~/ecalj/trash` へ（`47ce04cc5`）
- ecalj・ecaljdoc の控えのクローンを外付けドライブ `/media/takao/TAKAOMINI`（`ecalj_20261001`・`ecaljdoc_20261001`）に
- Materials Project の API キーを追跡から外した。`<ecalj>/MaterialProject.key`（無視される）と見本 `.example`、`pylib/mpkey.py`（`cb9b2d7b7`）
- GW1500: 物質ごとの注記の表 [`ecalj_auto/gw1500_notes_20261001.tsv`](../../ecalj_auto/gw1500_notes_20261001.tsv)、分類と条件の違い（`GW1500_status.md` §5、`b0fdfda19`）。
  5 月に収束した GOOD はやり直さない（run3 のキューを NOTCONV だけに）。Rb8（mp-1179832）は不正な構造と分かり止めた（`221e0e307`）
- `gwscconv`: ギャップの無い反復（金属）は固有値の変化で収束を判定（`1730f7e7e`）
- kt1・kr7（nvfortran GPU）・mic（ifx）で試験の組 inputs・mlo・afsym・install・samples がすべて PASS（HEAD `b695fa65c`）

### 2026-09-30

- `fix-idu10` をマージ（`f14270e83`）: `idu = 10 + mode` が `sigm` の無いとき二重計数なしの LDA+U になっていた（2023-09 から）。
  MLOsamples の SmP・GdCo5 の LDA+U を作り直した（`6e5817ffe`、`23429f0a3`）
- GW1500: FAILED の 147 物質を fp32 で回し直し、144 が収束（`70f79b23b`）
- Samples/Legacy を片付けて組み直し（AtomDimer/N2 など）、ecalj と ecaljdoc を clone し直した
