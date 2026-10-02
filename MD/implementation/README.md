# MD/implementation — GW の実装の覚え書き

GW の主な段（χ0 と W、Σ）の実装の覚え書き。もとは `ecaljdoc/implementation/`（2026-10-02 に MD/ へ移した）。
ファイルには 2026-05 の TOML 化より前に書いた部分がある。プログラム名とファイル名は、いまのコードで確かめる
（W と Σ はいま `hgw` の一つのプロセスで計算する。[ecaljclaude.md](../ecaljclaude.md) の hgw_combined の節）。

| ファイル | 中身 |
|---|---|
| [`hx0fp.md`](hx0fp.md) | RPA の応答関数 χ0 と W の実装（`hx0fp0`・`hrcxq` の頃） |
| [`hsfp0.md`](hsfp0.md) | 自己エネルギー Σ の実装（`hsfp0_sc`） |
| [`issue.md`](issue.md) | 検討事項（χ0 の周波数の平滑化など） |

関係するもの:

- コードの大局: [module_map.md](../module_map.md)
- 有限温度と均し方（`t_*`）: [kBT/README.md](../kBT/README.md)、ecaljdoc [manual/kBT.md](../../ecaljdoc/manual/kBT.md)
- 入力: ecaljdoc [manual/gwinput.md](../../ecaljdoc/manual/gwinput.md)、[manual/gwsc.md](../../ecaljdoc/manual/gwsc.md)
