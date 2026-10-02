# 更新履歴（要約）

ecalj の新しい機能と、結果が変わる修正の要約。新しい順。詳しい記録（全部の変更）は ecalj の
[Changes.txt](https://github.com/tkotani/ecalj/blob/main/Changes.txt)。各項目の説明はリンク先のページ。

## 最近の大きな話題（2026-09〜10）

### 1. 均し方のキー `t_*` と、LiTi₂O₄ の反復の結果（tf32・fp32・fp64）

- **温度（均し方）を `[gw]` の `t_*` のキーで指定する**（単位は K）。χ0 側は `t_tetrakbt`: 0 は T = 0 のテトラヘドロン法、正の値は有限温度のテトラヘドロン法、
  負の値（−T）は Im χ0 を温度 T の Fermi–Dirac と同じ幅の Gaussian で均す（以前の `SmearX0` に代わる。`SmearX0` が残っていると止まり、
  `ctrlg_update.py` が書き換える）。Σ 側は `t_sigmaw`（既定 1000 K、Σ の中間準位を均す Fermi–Dirac の核）。金属で反復や k メッシュによって
  バンドが荒れる原因（小さい q の W のプラズモン極を Σc の実軸の極の項が拾う）は `wcsmear`（既定）で均す → [kBT](./kBT) §0・§2・§3・§3.5
- **LiTi₂O₄（金属スピネル、14 原子）の MLO-QSGW を LDA から 40 反復**（6³ と 9³、`gwsc --mlo`）。GPU の行列積の精度（`gwsc --prec=tf32|fp32|fp64`）を変えても、
  MLO バンドの差は数 meV 以下（6³ の 10 反復で tf32 と fp32 が t2g で 0.8 meV、1 反復目の fp64 と fp32 が 1 meV）で、反復やメッシュによる差（数十 meV）より
  1 桁小さい。速い tf32 で回してよい（1 反復: 6³ で tf32 3.7 分、fp64 17 分。kt1、RTX 5090 × 2）。反復を 10 回で止めると収束の途中で、30〜40 反復が要る
  → [LiTi₂O₄ のまとめ](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/LiTi2O4/README.md)（表 3〜5、図 1・2）、[MLO-gwsc](./mlo_gwsc)、[GPU version](./ecaljgpu)

### 2. MLO とマグノン

- **MLO は Löwdin で直交化した関数（射影 Wannier 関数）を標準にした**（2026-10-02）。軌道の名前と対称性を保ち、模型は直交した基底の $H(\mathbf R)$ だけ。
  $U$・$J$・cRPA（`job_mloW`）とマグノン（`job_mlo_magnon`）もこの基底で、Wannier 関数（最大局在化）の経路は外した → [MLO](./mlo) §6
- **マグノン**（bcc Fe、FeCo、Ni）: MLO の窓によらなくなり、Wannier 関数による以前の計算に近づいた。bcc Fe では q ≥ 0.4 でほぼ重なり、q ≤ 0.3 では
  MLO が 2 割ほど高い。Goldstone の条件のための $W$ の倍率は Fe 1.24、FeCo 1.28、Ni 1.75
  → [Fe のマグノンのサンプル](https://github.com/tkotani/ecalj/blob/main/Samples/Magnon/Fe_mlo_magnon/README.md)（表 2・図 1）
- **cRPA の $U$**: Ni の d で MLO 2.90 eV、Wannier 3.78 eV（部分空間の取り方の違い。遮蔽をすべて入れた RPA の $U$ は 1.43 と 1.58 eV で近い）→ [MLO](./mlo) §6

## 2026-10-02

- **MLO は Löwdin で直交化した関数になった**。MLO の部分空間の射影 Wannier 関数で、軌道の名前（t₂g、e_g など）と対称性を保つ。
  模型は直交した基底の $H(\mathbf R)$ だけで持ち、スピン軌道、$U$・$J$・cRPA（`job_mloW`）、マグノン（`job_mlo_magnon`）、MLO-QSGW の $\Sigma$ もこの基底。
  メッシュ上のバンドは変わらず、メッシュの外の内挿は $E_F$ 近くで良くなる（63 物質で最大の誤差の中央値 0.019 → 0.011 eV）。
  マグノンは MLO の窓（`mlo_delta`、`mlo_w`）によらなくなり、Wannier 関数による以前の計算に近づいた。
  → [MLO](./mlo) §6「MLO の標準」、表 M8
  - **古い `HamRsMLO` は読まない**。`job_mlo` を回し直す。直交化しない以前の模型は `job_mlo <sname> --mlo_raw`（比べるときだけ）
  - 模型の検査（`mlo_bandcheck.py`）の判定 3 は「帯のとげ」（模型のバンドが 1 点で跳ぶ。一次従属の崩れ）になった → [MLO](./mlo) §9 式 (12)
- **対称性は spglib から**（`symgrp = "find"` のとき。spglib 2.6.0 を同梱、Python は要らない）。求めた操作は作業ディレクトリの
  `symmetry.<sname>.json` に書き、次からはそれを読む。超格子の純粋な並進、反強磁性（`[[site]]` の `af`）の磁気対称性も含む。
  `symgrp` に生成元を書いた入力は従来どおり → [lmf](./lmf) の SYMGRP
  - 操作の並びが前と違うので、並びをそのまま書き出す出力（BoltzTraP の `.struct` など）は並びが変わる
- **反強磁性の対称性は `symgrpaf = "find"`** と書けば、spglib が `af` の印から求める（生成元を書かなくてよい）。`symgrpaf` を使うときは `pwmode = 11`
  → [UsageDetailed](./UsageDetailed) の Antiferro symmetry
- **Wannier 関数（最大局在化）の経路、異常ホール伝導度（AHC）、`lmfham2` を外した**。模型は MLO だけ。cRPA は `job_mloW <sname> --crpa`、
  広がりは `mlo_spread.py`、マグノンは `job_mlo_magnon`。外す前の版は git のタグ `last-wannier` → [MLO](./mlo) §6
- MLO を実空間で規格化した（$U$・$J$ が長さ 1 の軌道の値になる。Fe の d の $W$ が 1〜6 % 動いた）→ [MLO](./mlo) §6 式 (7a)
- 直した誤り: `lmchk --getwsr`（MT 半径の見積もり）がスピン分極の入力で止まっていた

## 2026-10-01

- **MLO の模型の既定を整理した**（基準 1・2・3）。半内殻の局所軌道は帯の位置から自動で選ぶ（入力のキーは無い）。gwinit は 4f の原子の f と、
  `!` 付きの `mlo_lm2`（外せば基準 2、空隙の大きい MgTe・ZnTe・CdTe・AlN などが合う）を書く → [MLO](./mlo) §1・§9
- **模型の検査** `mlo_bandcheck.py`（`job_mlo`・`job_mlo_soc` が最後に回し、CHECK PASS/FAIL を出す）→ [MLO](./mlo) §9
- `job_mlo_soc` は、スピン軌道ありの DFT のバンドも描いて比べる → [MLO](./mlo) §4
- Materials Project の API キーは、ecalj の最上位の `MaterialProject.key`（git で無視）か環境変数 `MP_API_KEY` に置く

## 2026-09-30

- **結果が変わる修正**:
  - `pwmode = 11`（既定）での LDA+U の密度行列、スピン軌道の軌道モーメント、`--cls` が、別の k の次元で固有ベクトルを読んでいた（2026-03-30 から）。
    この間の `pwmode = 11` の LDA+U の SCF は誤り
  - `[[spec]]` の `idu = 10 + mode`（4f の `idu = 12` など）は、`sigm` が無いとき二重計数の項を引かない LDA+U になっていた（2023-09 から）
  - 反強磁性の対称性（`symgrpaf`）の QSGW が止まっていた → 両方のスピンを計算して動く
- `--ctrlg:<path>=<value>` で、ファイルに無いキーも書き足して効かせる（これまで既定値のまま走った）→ [TOML migration](./toml_migration)
- `gwscconv`: 金属・半金属は固有値の変化で収束を判定する → [gwsc](./gwsc)
- サンプルを組み直した（SOC の磁気異方性、LDA+U、有効質量、構造緩和、DOS、IIR、箱の中の分子、電子温度のスキャンなど）→ [Samples](./samples)

## 2026-09 以前（大きな話題）

| 時期 | 話題 | 説明 |
| --- | --- | --- |
| 2026-09 | MLO-QSGW（自己エネルギーを MLO 表現で内挿する QSGW、`gwsc --mlo`） | [MLO-gwsc](./mlo_gwsc) |
| 2026-09 | 有限温度（χ0 の `t_tetrakbt`、Σ の準位の幅 `t_sigmaw`、`wcsmear`） | [kBT](./kBT) |
| 2026-09 | GW の GPU 高速化（精度の切り替え `--prec=tf32\|fp32\|fp64`、方法は表で自動） | [GPU version](./ecaljgpu) |
| 2026-09 | MLO の自動の窓（`mlo_method = 4`）、局所軌道の扱い | [MLO](./mlo) |
| 2026-06 | 有限温度の四面体法、`gw_lmfh` の GPU | [kBT](./kBT)、[gwsc](./gwsc) |
| 2026-05 | 入力は `ctrlg.<sname>.toml` 一本に（TOML）、上書きは `--ctrlg:` | [TOML migration](./toml_migration) |
