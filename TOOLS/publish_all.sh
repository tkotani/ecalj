#!/bin/bash
# Publish everything in order (2026-10-03, for the maintainers):
#   1. ecalj: git push dev main (and rel main with --rel)
#   2. ecaljdoc: TOOLS/publish_ecaljdoc.sh (github.com/ecalj/ecaljdoc, the site https://ecalj.github.io/ecaljdoc/)
#   3. the GW1500 database: TOOLS/publish_gw1500db.sh (github.com/tkotani/DOSnpSupplement, QSGW80_2026/)
# The order matters: the site and the database link to files of the ecalj repository on GitHub.
#
#   TOOLS/publish_all.sh                 show what would be published (commits of ecalj ahead of dev/rel; ecaljdoc and the database
#                                        are copied into their publishing clones and committed there, not pushed)
#   TOOLS/publish_all.sh --push          push dev, publish ecaljdoc, push the database
#   TOOLS/publish_all.sh --push --rel    the same, and push rel main as well (after the tests on three machines, MD/ForDevelopers.md §1)
set -eu
TOP=$(git -C "$(dirname "$(readlink -f "$0")")" rev-parse --show-toplevel)
PUSH=no; REL=no
for a in "$@"; do case $a in --push) PUSH=yes;; --rel) REL=yes;; *) echo "unknown option $a"; exit 1;; esac; done
cd "$TOP"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "publish_all: the work tree of ecalj has uncommitted changes; commit them first"; git status --short | head; exit 1
fi
git fetch -q dev; [ $REL = yes ] && git fetch -q rel
echo "== 1. ecalj $(git rev-parse --short=9 HEAD): $(git rev-list --count dev/main..HEAD) commits ahead of dev/main$( [ $REL = yes ] && echo ", $(git rev-list --count rel/main..HEAD) ahead of rel/main")"
if [ $PUSH = yes ]; then
  git push dev main
  [ $REL = yes ] && git push rel main
fi
echo "== 2. ecaljdoc"
if [ $PUSH = yes ]; then "$TOP/TOOLS/publish_ecaljdoc.sh" --push; else "$TOP/TOOLS/publish_ecaljdoc.sh"; fi
echo "== 3. GW1500 database"
if [ $PUSH = yes ]; then "$TOP/TOOLS/publish_gw1500db.sh" --push; else "$TOP/TOOLS/publish_gw1500db.sh"; fi
echo "publish_all: done ($( [ $PUSH = yes ] && echo pushed || echo 'nothing pushed; add --push'))"
