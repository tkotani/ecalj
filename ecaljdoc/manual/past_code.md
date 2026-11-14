# 外したコード（過去のコード） {#past-code}

過去に使っていて、いまは使われなくなったコードの一覧です。リポジトリからは外しましたが、git の履歴に残っているので、
下の表の「外す前のコミット」から取り出せます。また使うかもしれないもの（試作など）も、ここから探してください。

```bash
git show <コミット>:<パス>                 # 1 つのファイルを表示する
git checkout <コミット> -- <パス>          # 1 つのファイルを作業ツリーに戻す
git show <コミット> --stat                 # そのコミットで何が変わったか
```

取り出しは、ecalj をクローンした手元で行います。開発者向けのもっと詳しい記録（何をどう判断したか）は、リポジトリの `MD/past_log.md` の表 1 にあります（このサイトには無い）。

*表 1*. 外したコード（新しいものが上）

| 外した日 | 何か | 元の場所 | 外す前のコミット |
| --- | --- | --- | --- |
| 2026-11-14 | どこからも使われなくなった module 5 つ: `m_bz_integ`（四面体の頂点の平均による BZ 積分）、`m_pbindex`（積基底の添字の表）、`m_readqgcou`（QGcou の読み込み）、`m_HamRsMLO`（古い `HamRsMLO` の読み込み）、`m_qdwh_eigensolver`（GPU の固有値ソルバの試作。QDWH 極分解による一般化固有値問題。2026-03。また使うかもしれない） | `SRC/subroutines/m_bz_intg.f90`、`m_pbindex.f90`、`m_readQGcou.f90`、`m_qdwh_eigensolver.f90`、`m_HamPMT.f90` の末尾の `module m_HamRsMLO` | `31d8516dc` |
| 2026-11-14 | GWinput を読む入力ルーチン `m_keyvalue`（`getkeyvalue`）。2026-05 の TOML 化で使われなくなった | `SRC/subroutines/keyvalue.f90` | `31d8516dc` |
| 2026-10-07 | 2026-05 より前の入力の流れのスクリプト 3 本 | `SRC/exec/ctrlgenM1.py`、`auto_creplot.py`、`auto_job_mp.py` | `72967dec5` |
| 2026-10-02 | Wannier 関数（最大局在化）の経路: `hmaxloc`、`hpsig`、`wanplot`、`hsocmat`、Wannier のマグノン `hmagnon`、スクリプト `genMLWF`、`job_magnon` など。MLO に置き換えた | `SRC/main/`、`SRC/subroutines/`（`main_wannier_*`、`m_wannier_maxloc` など）、`SRC/exec/` | タグ `last-wannier`（`dbcd6e51d`）。そこで `Samples/WannierVsMLO/run.sh` が MLO との比較を再現する |
| 2026-10-02 | AHC（異常ホール伝導度）: `hahc`、`x0kf_ahc`、`job_AHC`、`hx0ahc.py` | `SRC/main/`、`SRC/subroutines/`、`SRC/exec/` | タグ `last-wannier` |
| 2026-10-02 | lmfham2（反復で作る古い MLO）: `lmfham2`、`genMLO` | `SRC/main/`、`SRC/subroutines/`、`SRC/exec/` | タグ `last-wannier` |
| 2026-10-02 | 対称性のルーチン一式 `m_mksym_util.f90`（中身は `m_symop_util`、`m_symderive`、`m_symfind`、`m_mksym` に移した） | `SRC/subroutines/m_mksym_util.f90` | `b889a00b8` |
| 2026-10-01 | 退役したスクリプト | `SRC/exec_legacy/` | `a1c30d1a1` |
| 2026-10-01 | ジョブの自動実行の試作 | `jobauto/` | `a1c30d1a1` |
| 2026-10-01 | TOOLS の古い道具（約 60 項目） | `TOOLS/` | `c6a1b6b7f` |
