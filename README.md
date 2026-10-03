# ecalj

First-principles electronic-structure package (PMT = APW + MTO basis; LDA/GGA and QSGW).

- Manual: https://ecalj.github.io/ecaljdoc/
- With AI (New!): https://ecalj.github.io/ecaljdoc/manual/getstartedAI — give this page to an AI assistant to be guided through ecalj
- Start here: [MD/README.md](MD/README.md)
- Changes: [Changes.md](Changes.md)

## For maintainers: publishing (push)

Three repositories, in this order (the site and the database link to files of ecalj on GitHub):

1. ecalj: `git push dev main`, and `git push rel main` after the tests on three machines (`MD/ForDevelopers.md` §1)
2. the manual: `TOOLS/publish_ecaljdoc.sh --push` → github.com/ecalj/ecaljdoc, the site https://ecalj.github.io/ecaljdoc/ (`MD/ecaljdoc_publish.md`)
3. the GW1500 database: `TOOLS/publish_gw1500db.sh --push` → github.com/tkotani/DOSnpSupplement, directory `QSGW80_2026/`
   (built by `ecalj_auto/gw1500db_build.py`; the 2025 supplement of arXiv:2507.19189 there is kept)

All at once: `TOOLS/publish_all.sh` shows what would be published (ecaljdoc and the database are committed in their publishing
clones, nothing is pushed); `TOOLS/publish_all.sh --push` pushes dev and publishes 2 and 3; add `--rel` to push rel as well.
