# WannierVsMLO: MLO と Wannier 関数（最大局在化）の比較

2026-10-02 に ecalj から Wannier 関数の経路（`genMLWFx`、`hmaxloc`、`hpsig_MPI` など）を外す前に、MLO で置き換えた cRPA と広がりを比べた記録。
**この比較は、git のタグ `last-wannier` のコミットでしか回らない**（その後のコミットには Wannier の経路が無い）。

```bash
git checkout last-wannier          # ビルドし直す（InstallAll.py）
cd Samples/WannierVsMLO && ./run.sh -np 8
```

- 入力は `TestInstall/ni_crpa`・`srvo3_crpa`（当時の Wannier の試験）に `[mlo] mlo_nkabc` を足し、SrVO₃ は GW の k メッシュを 4³ にしたもの。
  模型の部分空間は `mlo_lm`（Ni d 5 本、V t₂g 3 本）で、Wannier も同じ `mlo_lm` を読む
- MLO: `job_band` → `job_mlo` → `job_mloW --crpa` → `mlo_spread.py`。Wannier: `job_band` → `genMLWFx`（`bnds.<sname>` が要る）
- `ujk.py`: `Coulomb_v.UP`、`Screening_W-v.UP`、`Screening_W-v_crpa.UP` から ω = 0、R = 0 の U = (ii|ii)、U′ = (ii|jj)、J = (ij|ji) の平均
- `wan_spread.py`: Wannier 関数の広がりを MLO と同じ式（ecaljdoc mlo の式 (7e)）で。`UUU`（バンドの間の ⟨u_k|u_k+b⟩）と `MLWU`（ゲージ行列）から。
  Marzari–Vanderbilt の形も出し、`hmaxloc` の値と一致することを確かめた

結果（2026-10-02、t14 gfortran、-np 4）は ecaljdoc [mlo](https://ecalj.github.io/ecaljdoc/manual/mlo) §6 の表 M7。
SrVO₃（4³）は cRPA の U が MLO 3.125 eV、Wannier 3.149 eV。Ni d は MLO 2.84 eV、Wannier 3.78 eV（MLO の d が広く、除く遮蔽が少ない）。
細かい数値は研究ログ `MD/research_log.md` 2026-10-02 00:32。
