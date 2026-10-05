#!/bin/bash
# Publish the GW1500 QSGW80 database to github.com/tkotani/DOSnpSupplement, in its subdirectory QSGW80_2026/ (2026-10-03).
# The 2025 content of that repository (the supplement of arXiv:2507.19189) is left as it is; its top README gets one link line.
#
#   TOOLS/publish_gw1500db.sh          copy the database into the publishing clone and commit there; show the change
#   TOOLS/publish_gw1500db.sh --push   ... and push it (only when the maintainer says so)
#
# Env: GW1500DB (default /media/takao/TAKAOMINI/gw1500db: made by ecalj_auto/gw1500db_build.py),
#      DOSNP_PUBLISH_DIR (default ~/work/DOSnpSupplement_publish, a clone, made if missing),
#      DOSNP_REMOTE (default git@github.com:tkotani/DOSnpSupplement.git).
# Published: README.md, table.md, bands_*.md, fig/, gw1500db.tsv, history.tsv, summary_counts.json. Not published: npz/ and logs/ (raw data,
# kept by the maintainers), the side files of the builder.
set -eu
TOP=$(git -C "$(dirname "$(readlink -f "$0")")" rev-parse --show-toplevel)
DB=${GW1500DB:-/media/takao/TAKAOMINI/gw1500db}
PUB=${DOSNP_PUBLISH_DIR:-$HOME/work/DOSnpSupplement_publish}
REMOTE=${DOSNP_REMOTE:-git@github.com:tkotani/DOSnpSupplement.git}
SUB=QSGW80_2026
PUSH=no; [ "${1:-}" = "--push" ] && PUSH=yes
[ -f "$DB/README.md" ] && [ -d "$DB/fig" ] || { echo "publish_gw1500db: no database at $DB (run ecalj_auto/gw1500db_build.py <db> --figs)"; exit 1; }
REV=$(git -C "$TOP" rev-parse --short=9 HEAD)

[ -d "$PUB/.git" ] || git clone -q "$REMOTE" "$PUB"
git -C "$PUB" fetch -q origin
git -C "$PUB" checkout -q main
git -C "$PUB" reset -q --hard origin/main

mkdir -p "$PUB/$SUB"
rsync -a --delete --include='README.md' --include='table.md' --include='bands_*.md' --include='gw1500db.tsv' \
  --include='summary_counts.json' --include='history.tsv' --include='fig/' --include='fig/*.png' --exclude='*' "$DB/" "$PUB/$SUB/"

# one link line in the top README (once)
LINK="**New (2026): [QSGW80 band gaps, bands and DOS of 1546 materials, iterated to convergence]($SUB/README.md)** — the GW1500 database of ecalj; the 2025 tables below are kept as the supplement of arXiv:2507.19189."
grep -q "($SUB/README.md)" "$PUB/README.md" || { printf '%s\n\n' "$LINK" | cat - "$PUB/README.md" > "$PUB/README.md.new" && mv "$PUB/README.md.new" "$PUB/README.md"; }

git -C "$PUB" add -A
if git -C "$PUB" diff --cached --quiet; then
  echo "publish_gw1500db: nothing to publish (the database in $DB is the same as $REMOTE main)"; exit 0
fi
when=$(grep -m1 -o 'Made [0-9-]* [0-9:]*' "$DB/README.md" | sed 's/Made //')
git -C "$PUB" commit -q -m "QSGW80_2026: GW1500 database of $when (ecalj $REV)"
echo "publish_gw1500db: committed in $PUB:"
git -C "$PUB" show --stat --format='%h %s' HEAD | tail -5
if [ "$PUSH" = yes ]; then
  git -C "$PUB" push origin main
  echo "publish_gw1500db: pushed to $REMOTE"
else
  echo "publish_gw1500db: not pushed. To publish: TOOLS/publish_gw1500db.sh --push   (or git -C $PUB push origin main)"
fi
