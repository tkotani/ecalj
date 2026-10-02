# Samples/AtomDimer — 箱の中の分子と原子

> 説明の本体: ecaljdoc の [manual/UsageDetailed.md](../../ecaljdoc/manual/UsageDetailed.md)「For molecules, we may use --systype=molecule for ctrlgenToml.py.」（サイト https://ecalj.github.io/ecaljdoc/manual/UsageDetailed）。サンプルの一覧は [manual/samples.md](../../ecaljdoc/manual/samples.md)。

| ディレクトリ | 内容 | 試験 |
| --- | --- | --- |
| [`N2/`](N2/README.md) | N₂ 分子と N 原子を 15 Å の箱で（PBE、スピン分極、固定磁気モーメント）。結合長 3 点から平衡結合長と結合エネルギー | `testecalj N2 -np 8`（7 分） |
| `elements_2012.txt` | 2011〜2012 年に周期表の 36 元素の同核二原子分子と原子を回したときの設定の記録（基底、初期の結合長、`fsmom`、箱） | — |

ほかの元素の分子は、`N2/` の入力の元素・`fsmom`・`mmom` を変えて作る（`N2/README.md`）。
