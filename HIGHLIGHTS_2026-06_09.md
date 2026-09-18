# ecalj 2026-06 〜 2026-09 の更新ハイライト

2026-06-03 以降、main に入った 140 コミットの要点。詳細は `Changes.txt`（日付順）と
[ecaljdoc](https://ecalj.github.io/ecaljdoc/)。

## 1. 入力は `ctrlg.<sname>.toml` 一つになった（破壊的変更）

- lmf も GW も **`ctrlg.<sname>.toml` だけを読む**。`PB.<sname>.toml`（product basis の表）、
  `esm_input.dat`（ESM）、`GWinput.toml` は廃止。残っていると abort して変換コマンドを案内する。
- セクション順は `[gw]` `[mlo]` `[blocks]` `[product_basis]`（末尾に nlx / valence / core の表）。
  `Worb` → `[mlo] mlo_lm`、`QforEPS` / `QforGW` → `[gw]`。
- 体裁: 空行は `# === X ===` の大見出しの前だけ、セクション内は `# ----` の罫線で区切る。
  `ctrlgenToml.py` / `gwinit` / `Legacy2toml.py` が同じ形を出す。`[esm]` は見え消しのテンプレが入る。
- **古い作業ディレクトリを使うとき**:
  ```
  ctrlg_absorb.py <sname>     # PB.<sname>.toml / esm_input.dat を ctrlg に取り込む（原本は .bk）
  Legacy2toml.py <sname>      # 旧 ctrl.<sname> + GWinput から
  ```
- 新規は `ctrls.<sname>`（構造だけ）→ `ctrlgenToml.py <sname>` → `lmfa` / `lmf` / `gwsc`。
  [manual/lmf](https://ecalj.github.io/ecaljdoc/manual/lmf)、[toml_migration](https://ecalj.github.io/ecaljdoc/manual/toml_migration)。

## 2. MLO（局在軌道）が実用形に

- 既定を `mlo_method = 4` に。パラメタは `mlo_delta` / `mlo_w`（eV, 既定 2.0）の 2 つだけ。
  [manual/mlo](https://ecalj.github.io/ecaljdoc/manual/mlo) に理論・損失関数・走査結果。
- `mlo_lm` で lm チャネルを原子ごとに指定（`5 6 8` で t2g だけ、など）。
  旧 Worb は殻全体しか選べていなかったバグを修正。
- `Samples/MLOsamples/` 25 系（Si, GaAs, C, Cu, Fe, NiO, SrTiO3, Al2O3:Cr, GdCo5, RuO2, SmP, FeMgO …、
  Materials Project から既定のまま回した Ag, Al, NaCl, SiC, CdTe, ZnO, TiO2）すべて method 4 で
  `testecalj` の回帰チェック。全系のフィッティング図は [manual/mlo](https://ecalj.github.io/ecaljdoc/manual/mlo)。
- FeMgO スラブは空格子球 + ESM（`[esm]` セクション）の 76 軌道モデルを標準に。

## 3. 有限温度 QSGW（試験的）

- `tetrakbt` / `t_tetrakbt`（χ₀ 側の Fermi–Dirac 占有と有限温度 E_F）と
  `t_sigmakbt`（Σ 側）。**既定は off**。`Samples/kBT/`（LiTi2O4 6³/9³, 1000–3000 K, Fe の対照実験）。
- [manual/kBT](https://ecalj.github.io/ecaljdoc/manual/kBT): 手法、実装レビュー（§7）、
  虚時間を使わない定式化（§8）、**残る課題（§9）**。300 K でも Σ 側を温めると金属で
  Σ−v_xc が 0.1 eV 動くなど、使うには §9 を読んでから。

## 4. 一発 GW と GPU

- `gw_lmfh <sname> -np N --gpu --mp` で G₀W₀ 型の QP エネルギー。`hsfp0_gpu`（W 縮約の CUDA 化）。
- 診断: `--dumpW`（W を保存）、`--WVR2ptRaxis`。`FiniteT_and_QPE_HOWTO.md`。

## 5. 直したバグ

- `lmf --quit=band` が LDA+U の `dmats.<sname>` を NaN で壊していた。
- `zhev`: LAPACK の ier を先に見て原因を名指しする。`zhev_tk2` の nev を n に丸める。
- PROCAR の k 点順序（np ≥ 11 で rank 接尾辞の数値ソート）。
- `m_sxcf_sc` の実軸極ビニングの範囲外アクセス。`readbandedge` の黙った fallback。
- heftet が絶縁体で `EFERMI_kbt` を書かず、有限温度モードの絶縁体が落ちていた。
- MLO の `nskip`（射影子から外す半芯状態の数）を k ごとに決めていたため、Cu の d 模型のように
  最下位が s 帯になる k とならない k で射影子が入れ替わり、バンドに折れが出ていた。
  全 k での最小値に固定し、外した帯と残した帯の間にギャップが無ければ止まる。
- `lmf` が前の run の混合履歴 `__mixm.<sname>` を引き継いでいた。非磁性に収束した run の履歴を
  継ぐと Broyden の 1 歩目で非磁性解へ落ちる（Fe 2.13 → 0.02 μB）。起動時に捨てるようにした
  （`--keepmixm` で従来どおり）。
- Zn 3d のような浅い半芯 LO（`pz`）を `mlo_lm` の模型関数に使えず（常に EH 関数だった）、
  ZnO の MLO が 476 meV 外れていた。LO の帯が E_F−10 eV より上なら LO を使う。

## 6. サンプルとテスト

- `testecalj --all` は 25 ターゲット（fe, gdn を追加）。`Samples/Legacy/` 32 ディレクトリ 53 ファイルを
  TOML 化、AFsymmetry の NiO/NiSe を現行文法に。`job_magnon --sp1/--sp2`。
- `Samples/GetStarted/GaAs/` が構造 → QSGW バンドの最短経路。
- 検証: ローカル gfortran-14 と kt1 nvfortran GPU+MP（clean clone）で
  `testecalj --all` 25/25、MLO 18/18、kBT Fe、GaAs 一気通貫。ifx（ucgw）は未実施。
