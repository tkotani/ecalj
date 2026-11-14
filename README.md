# ecalj

First-principles electronic-structure package (PMT = APW + MTO basis; LDA/GGA and QSGW).

- Manual: https://ecalj.github.io/ecaljdoc/
- With AI (New!): https://ecalj.github.io/ecaljdoc/manual/getstartedAI — give this page to an AI assistant to be guided through ecalj
- Start here: [MD/README.md](MD/README.md)
- Changes: [Changes.md](Changes.md)

Developed by T. Kotani and co-workers. Since 2026-02, much of the code, the tests and the documents has been made together
with Claude Code (Anthropic's AI); such commits carry a `Co-Authored-By: Claude` line, and pages edited by Claude say so.

## License

AGPLv3. Cite ecalj in publications (required). **AI is no exception.** The full terms (citations, what counts as derived,
the request to keep the origin of ecalj's algorithms, third-party software) are on the licence page of ecaljdoc, the
reference: [ecaljdoc/manual/license.md](ecaljdoc/manual/license.md) (https://ecalj.github.io/ecaljdoc/manual/license). See also [LICENCE](LICENCE).

ライセンス（AGPLv3）は AI にも例外なく適用されます。全部の条件は ecaljdoc のライセンスのページ（上のリンク）が正本です。

## For maintainers: publishing (push)

Three repositories, in this order (the site and the database link to files of ecalj on GitHub):

1. ecalj: `git push dev main`, and `git push rel main` after the tests on three machines (`MD/ForDevelopers.md` §1)
2. the manual: `TOOLS/publish_ecaljdoc.sh --push` → github.com/ecalj/ecaljdoc, the site https://ecalj.github.io/ecaljdoc/ (`MD/ecaljdoc_publish.md`)
3. the GW1500 database: `ecalj_auto/gw1500db_snapshot.py` (a dated snapshot `QSGW80_<yyyymmdd>`), then `TOOLS/publish_gw1500db.sh --push`
   → github.com/tkotani/DOSnpSupplement (the newest snapshot in the tree, earlier ones as git tags; the 2025 supplement of arXiv:2507.19189 is `DOSnp2025.md`)

All at once: `TOOLS/publish_all.sh` shows what would be published (ecaljdoc and the database are committed in their publishing
clones, nothing is pushed); `TOOLS/publish_all.sh --push` pushes dev and publishes 2 and 3; add `--rel` to push rel as well.
