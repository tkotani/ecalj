# gadget — 標準の流れに入れない、試作の道具

`InstallAll.py` は bindir に入れない。使うときはここから直接回す（`python3 <ecalj>/TOOLS/gadget/<名前>`）。
メインの機能ではないが、捨てると作り直すのに手間がかかり、何かに役に立ちうるものを置く（user 2026-10-02）。

## MLO の部分空間の中での最大局在化（Marzari–Vanderbilt）

2026-10-02 に試作。ecalj の MLO の標準は Löwdin で直交化した関数（ecaljdoc mlo §6）で、MV はメインにしない（user 2026-10-02）。

| 道具 | 何をするか |
| --- | --- |
| `mlo_maxloc.py <sname> [--sym] [--orb 5,6,7,8,9] [--save-u U.npz]` | `mlo_spread.py`（`SRC/exec`）の M(k,b) から、MLO・Löwdin・Löwdin + MV の広がり Ω と Ω_I を出す。`--sym` は Sakuma（PRB 87, 235109）の点群の拘束 |
| `mlo_cmlo_transform.py mv U_isp1.npz [U_isp2.npz]` | `__cmlo.data` を MV の基底に書き換える（`job_mloW`・`job_mlo_magnon` の W や K をその基底で取るため）。`lowdin` は今は恒等変換 |

役に立ちうる所:
- 共有結合の系（Si・C の sp³ など）で、原子の上でなく結合の上の局在関数が要るとき（拘束なしの MV は結合の方向にずれた混成軌道を作る）
- Wannier90 などの最大局在 Wannier 関数と、同じ部分空間で比べるとき（U、J、マグノン）
- Löwdin の関数がどれだけ最大局在に近いかの目安（Fe の t₂g で MV がさらに 7 % 縮めるだけ、Ni の d では Ω の 97 % が Ω_I で効かない。ecalj `MD/research_kotani_log.md` 2026-10-02 06:20）

制限: `--sym` は 1 原子・symmorphic な結晶だけ。`mlo_cmlo_transform.py mv` はメッシュの外の q（GW の q0、マグノンのずらしたメッシュ）では一番近いメッシュの点の U を使う。
`mlo_maxloc.py` は `SRC/exec` の `mlo_spread.py`・`symfind.py` を読み、`--sym` では PATH の上の `lmf` の場所（bindir）の `libecaljF.so` の `rotdlmm` を呼ぶ。
