#!/bin/bash
# Publish a dated snapshot of the GW1500 QSGW80 database to github.com/tkotani/DOSnpSupplement (2026-10-03; dated snapshots and one
# cover page 2026-10-07, user: "データベースごとに表紙を用意（メインの表紙を一枚のものにしてそこからリンクする）").
#
#   TOOLS/publish_gw1500db.sh [<snapshot dir>]          copy the snapshot into the publishing clone, write the cover, commit; show the change
#   TOOLS/publish_gw1500db.sh [<snapshot dir>] --push   ... and push it (only when the maintainer says so)
#
# <snapshot dir>: made by ecalj_auto/gw1500db_snapshot.py, named QSGW80_<yyyymmdd> (default: the newest in
#   /media/takao/TAKAOMINI/gw1500db_snapshots). It is published as the subdirectory of the same name; earlier snapshots
#   (QSGW80_2026 of 2026-10-03, ...) stay as they are.
# The cover README.md lists every QSGW80_* snapshot and the 2025 database. The 2025 page (the supplement of arXiv:2507.19189) moves
# from README.md to DOSnp2025.md once, unchanged; its files (GW/, LDA/, bandpng.md, ...) stay where they are.
# Env: DOSNP_PUBLISH_DIR (default ~/work/DOSnpSupplement_publish, a clone, made if missing),
#      DOSNP_REMOTE (default git@github.com:tkotani/DOSnpSupplement.git), SNAPDIR (default /media/takao/TAKAOMINI/gw1500db_snapshots).
set -eu
TOP=$(git -C "$(dirname "$(readlink -f "$0")")" rev-parse --show-toplevel)
PUB=${DOSNP_PUBLISH_DIR:-$HOME/work/DOSnpSupplement_publish}
REMOTE=${DOSNP_REMOTE:-git@github.com:tkotani/DOSnpSupplement.git}
SNAPS=${SNAPDIR:-/media/takao/TAKAOMINI/gw1500db_snapshots}
PUSH=no; SNAP=
for a in "$@"; do case "$a" in --push) PUSH=yes;; *) SNAP=$a;; esac; done
[ -n "$SNAP" ] || SNAP=$(ls -d "$SNAPS"/QSGW80_* 2>/dev/null | sort | tail -1)
[ -f "$SNAP/README.md" ] && [ -d "$SNAP/fig" ] || { echo "publish_gw1500db: no snapshot at '$SNAP' (ecalj_auto/gw1500db_snapshot.py)"; exit 1; }
SUB=$(basename "$SNAP")
REV=$(git -C "$TOP" rev-parse --short=9 HEAD)

[ -d "$PUB/.git" ] || git clone -q "$REMOTE" "$PUB"
git -C "$PUB" fetch -q origin
git -C "$PUB" checkout -q main
git -C "$PUB" reset -q --hard origin/main

mkdir -p "$PUB/$SUB"
rsync -a --delete "$SNAP/" "$PUB/$SUB/"

# the 2025 page out of README.md (once): the old link line of 2026-10-03 is dropped, the rest is kept as it was
if [ ! -f "$PUB/DOSnp2025.md" ]; then
  grep -v '^\*\*New (2026): \[QSGW80' "$PUB/README.md" | sed '1{/^$/d}' > "$PUB/DOSnp2025.md"
fi

# the cover
python3 - "$PUB" <<'EOF'
import sys, os, re, glob
pub = sys.argv[1]
rows = []
for d in sorted(glob.glob(f'{pub}/QSGW80_*/README.md'), reverse=True):
    sub = os.path.basename(os.path.dirname(d)); t = open(d).read()
    made = re.search(r'Made (\d{4}-\d\d-\d\d \d\d:\d\d)', t); made = made.group(1) if made else ''
    st = re.search(r'Models so far: (\d+) \(([^)]*)\)', t)
    mlo = f'MLO models: {st.group(2)}' if st else 'MLO models: —'
    files = 'tables, figures, history' + (', inputs (ctrlg), band data' if os.path.isdir(f'{pub}/{sub}/inputs') else '')
    rows.append(f'| [{sub}]({sub}/README.md) | {made} | QSGW80 (scaledsigma 0.8) iterated to convergence, LDA, MLO model; {mlo} | {files} |')
cover = f'''# DOSnpSupplement — band gaps of the GW1500 materials by ecalj

Databases of band gaps, band structures and DOS of about 1500 nonmagnetic materials of the Materials Project, computed with
[ecalj](https://github.com/tkotani/ecalj) (QSGW, LDA, MLO models). Each database has its own cover page with the conditions,
the version of ecalj of every value, and the checks. Newer databases do not replace older ones: values keep the conditions
and the code of their date.

| database | made | content | files |
| --- | --- | --- | --- |
''' + '\n'.join(rows) + '''
| [DOSnp2025](DOSnp2025.md) | 2025 | 1shot and 2ndshot QSGW and LDA of 1516 materials, the supplement of arXiv:2507.19189 | tables, band plots (GW/, LDA/) |

The crystal structures are those of the Materials Project (`ecalj_auto/INPUT/gw1500/POSCARALL` in ecalj).

Erratum (2026-10-07): in the figures of QSGW80_2026 the DOS panel is shifted by E_F (the gaps and the bands are right);
corrected from QSGW80_20261007 on.
'''
open(f'{pub}/README.md', 'w').write(cover)
EOF

git -C "$PUB" add -A
if git -C "$PUB" diff --cached --quiet; then
  echo "publish_gw1500db: nothing to publish ($SUB is the same as $REMOTE main)"; exit 0
fi
when=$(grep -m1 -o 'Made [0-9-]* [0-9:]*' "$SNAP/README.md" | sed 's/Made //')
git -C "$PUB" commit -q -m "$SUB: GW1500 database of $when (ecalj $REV); cover page"
echo "publish_gw1500db: committed in $PUB:"
git -C "$PUB" show --stat --format='%h %s' HEAD | tail -5
if [ "$PUSH" = yes ]; then
  git -C "$PUB" push origin main
  echo "publish_gw1500db: pushed to $REMOTE"
else
  echo "publish_gw1500db: not pushed. To publish: TOOLS/publish_gw1500db.sh $SNAP --push   (or git -C $PUB push origin main)"
fi
