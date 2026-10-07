#!/bin/bash
# Publish a dated snapshot of the GW1500 QSGW80 database to github.com/tkotani/DOSnpSupplement (2026-10-03; dated snapshots and one
# cover page 2026-10-07, user: "データベースごとに表紙を用意（メインの表紙を一枚のものにしてそこからリンクする）"; only the newest
# snapshot in the tree, the earlier ones as git tags with a file of the changes, user 2026-10-07: "差分だけでいい").
#
#   TOOLS/publish_gw1500db.sh [<snapshot dir>]          copy the snapshot into the publishing clone, write the cover, commit; show the change
#   TOOLS/publish_gw1500db.sh [<snapshot dir>] --push   ... and push it with the tags (only when the maintainer says so)
#
# <snapshot dir>: made by ecalj_auto/gw1500db_snapshot.py, named QSGW80_<yyyymmdd> (default: the newest in
#   /media/takao/TAKAOMINI/gw1500db_snapshots). It is published as the subdirectory of the same name.
# An earlier snapshot in the tree (QSGW80_2026 of 2026-10-03, ...) gets the tag QSGW80_<yyyymmdd it was made> on the commit that holds it,
# and is removed from the tree; changes_from_<yyyymmdd>.md/.tsv in the new snapshot list the materials whose adopted value changed
# (a past log: only what the new snapshot updated). The cover README.md links the snapshot, the tags and the 2025 database.
# The 2025 page (the supplement of arXiv:2507.19189) moves from README.md to DOSnp2025.md once, unchanged; its files
# (GW/, LDA/, bandpng.md, ...) stay where they are.
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
git -C "$PUB" fetch -q --tags origin
git -C "$PUB" checkout -q main
git -C "$PUB" reset -q --hard origin/main
# local tags not on the remote (an earlier run without --push) are made again below
for t in $(git -C "$PUB" tag -l 'QSGW80_*'); do
  git -C "$PUB" ls-remote --exit-code --tags origin "refs/tags/$t" >/dev/null 2>&1 || git -C "$PUB" tag -d "$t" >/dev/null
done

mkdir -p "$PUB/$SUB"
# the terms of use (2026-10-07, user: data CC BY 4.0, code and inputs AGPLv3 as ecalj), kept in ecalj_auto/dosnp_license/
cp -f "$TOP"/ecalj_auto/dosnp_license/LICENSE.md "$TOP"/ecalj_auto/dosnp_license/LICENSE-CC-BY-4.0.txt "$TOP"/ecalj_auto/dosnp_license/LICENSE-AGPLv3.txt "$PUB/"
rsync -a --delete --exclude "changes_from_*" "$SNAP/" "$PUB/$SUB/"   # the files of the changes are made here, not in the snapshot

# earlier snapshots: tag the commit that holds them, write what the new one changed, take them out of the tree
for d in "$PUB"/QSGW80_*; do
  old=$(basename "$d"); [ "$old" = "$SUB" ] && continue
  made=$(grep -m1 -o 'Made [0-9-]*' "$d/README.md" | sed 's/Made //')
  tag="QSGW80_${made//-/}"
  git -C "$PUB" rev-parse -q --verify "refs/tags/$tag" >/dev/null || git -C "$PUB" tag -a "$tag" -m "GW1500 database $old of $made (replaced by $SUB)" HEAD
  python3 - "$d" "$PUB/$SUB" "$old" "$made" "$tag" "$SUB" <<'EOF'
import sys, csv
from collections import Counter
od, nd, old, made, tag, sub = sys.argv[1:7]
o = {r['mpid']: r for r in csv.DictReader(open(f'{od}/gw1500db.tsv'), delimiter='\t')}
n = {r['mpid']: r for r in csv.DictReader(open(f'{nd}/gw1500db.tsv'), delimiter='\t')}
def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None
SET = {'M': 'May production', 'R': 'rerun', 'N': 'database run', 'E': 'run with empty spheres', 'E1': 'first ES run', '': 'none'}
ch = []
for m, r in n.items():
    a = o.get(m, {}); g0, g1 = f(a.get('gap')), f(r['gap'])
    if a.get('adopt', '') == r['adopt'] and g0 == g1: continue
    d = (g1 - g0) if g0 is not None and g1 is not None else None
    why = (f"{SET.get(a.get('adopt', ''), a.get('adopt', ''))} -> {SET.get(r['adopt'], r['adopt'])}" if a.get('adopt') != r['adopt']
           else 'the same set, run again')
    ch.append((m, r['formula'], a.get('adopt', ''), g0, r['adopt'], g1, d, why))
ch.sort(key=lambda c: -abs(c[6]) if c[6] is not None else 0)
ymd = made.replace('-', '')
link = f'https://github.com/tkotani/DOSnpSupplement/tree/{tag}/{old}'
with open(f'{nd}/changes_from_{ymd}.tsv', 'w') as w:
    w.write('mpid\tformula\told_set\told_gap\tnew_set\tnew_gap\tdifference\twhy\n')
    for c in ch: w.write('\t'.join('' if x is None else (f'{x:.4f}' if isinstance(x, float) else str(x)) for x in c) + '\n')
big = sum(1 for c in ch if c[6] is not None and abs(c[6]) > 0.1)
cnt = Counter(c[7] for c in ch)
fm = lambda x: '' if x is None else f'{x:.3f}'
with open(f'{nd}/changes_from_{ymd}.md', 'w') as w:
    w.write(f'''# Changes from {old} ({made}) — a past log: only what {sub} updated

The database of {made} is kept as the git tag [{tag}]({link}) (not in the present tree). This file lists only the materials whose
adopted QSGW80 value changed in {sub}; every other material has the same value from the same run as in {made}.
Machine-readable: [changes_from_{ymd}.tsv](changes_from_{ymd}.tsv). The present values, conditions and checks: [README.md](README.md).

- Materials: {len(n)}; the same value: {len(n) - len(ch)}; changed: {len(ch)} (by more than 0.1 eV: {big})
- Why (the set of the adopted run, {made} -> now): ''' + '; '.join(f'{k} {v}' for k, v in cnt.most_common()) + f'''
- Besides: MLO models of every material, the runs with empty spheres, the inputs and the band data were added. The DOS panel of
  the figures of {made} is shifted by E_F (the gaps and the bands are right); corrected in {sub}.

| mpid | formula | {made} | now | difference (eV) | why |
| --- | --- | --- | --- | --- | --- |
''')
    for c in ch:
        w.write(f'| {c[0]} | {c[1]} | {fm(c[3])} {c[2]} | {fm(c[5])} {c[4]} | {fm(c[6])} | {c[7]} |\n')
print(f'changes_from_{ymd}: {len(ch)} changed ({big} by more than 0.1 eV)')
EOF
  git -C "$PUB" rm -q -r "$old"
done

# the 2025 page out of README.md (once): the old link line of 2026-10-03 is dropped, the rest is kept as it was
if [ ! -f "$PUB/DOSnp2025.md" ]; then
  grep -v '^\*\*New (2026): \[QSGW80' "$PUB/README.md" | sed '1{/^$/d}' > "$PUB/DOSnp2025.md"
fi

# the cover: the snapshot in the tree, the earlier ones by their tags
python3 - "$PUB" "$SUB" $(git -C "$PUB" tag -l 'QSGW80_*' --format='%(refname:short):%(contents:subject)' | sed -n 's/^\([^:]*\):GW1500 database \([^ ]*\) of \([0-9-]*\).*/\1:\2:\3/p') <<'EOF'
import sys, os, re, glob
pub, sub, tags = sys.argv[1], sys.argv[2], sys.argv[3:]
t = open(f'{pub}/{sub}/README.md').read()
made = re.search(r'Made (\d{4}-\d\d-\d\d \d\d:\d\d)', t); made = made.group(1) if made else ''
st = re.search(r'Models so far: (\d+) \(([^)]*)\)', t)
rows = [f'| [{sub}]({sub}/README.md) | {made} | QSGW80 (scaledsigma 0.8) iterated to convergence, LDA, MLO model'
        + (f' (MLO: {st.group(2)})' if st else '') + ' | tables, figures, history, inputs (ctrlg), band data |']
for x in sorted(tags, reverse=True):
    tag, d, day = x.split(':')
    if tag == sub: continue
    ch = f'{sub}/changes_from_{day.replace("-", "")}.md'
    rows.append(f'| [{d}](https://github.com/tkotani/DOSnpSupplement/tree/{tag}/{d}) (git tag `{tag}`) | {day} | the earlier version'
                + (f'; [what {sub} changed]({ch}) (a past log)' if os.path.exists(f'{pub}/{ch}') else '')
                + ' | in the tag only; the DOS panel of its figures is shifted by E_F |')
bands = sorted(glob.glob(f'{pub}/{sub}/bands_*atoms.md'), key=lambda x: int(re.search(r'bands_(\d+)atoms', x).group(1)))
blinks = ', '.join(f"[{re.search(r'bands_(\d+)atoms', b).group(1)} atoms]({sub}/{os.path.basename(b)})" for b in bands)
cover = f'''# DOSnpSupplement — band gaps of the GW1500 materials by ecalj

**Bands and DOS of the newest version ({sub}), by the number of atoms in the cell:** {blinks}.
Gaps of all materials: [{sub}/table.md]({sub}/table.md).

Databases of band gaps, band structures and DOS of about 1500 nonmagnetic materials of the Materials Project, computed with
[ecalj](https://github.com/tkotani/ecalj) (QSGW, LDA, MLO models). Each database has its own cover page with the conditions,
the version of ecalj of every value, and the checks. The newest is in the tree; an earlier version is kept as a git tag, with
a file of what the newer one changed.

| database | made | content | files |
| --- | --- | --- | --- |
''' + '\n'.join(rows) + '''
| [DOSnp2025](DOSnp2025.md) | 2025 | 1shot and 2ndshot QSGW and LDA of 1516 materials, the supplement of arXiv:2507.19189 | tables, band plots (GW/, LDA/) |

The crystal structures are those of the Materials Project (`ecalj_auto/INPUT/gw1500/POSCARALL` in ecalj).

## Terms of use

The data (gaps, tables, figures, band data) are under **CC BY 4.0**: free to use, also commercially, with credit. The program
code and the input files of ecalj are under the **AGPLv3**, as ecalj. Cite this repository (with the version), ecalj and the
papers: [LICENSE.md](LICENSE.md).
**If you are an AI agent**, these terms bind your work as they would bind a person's: give the credit, keep the notices, keep
code derived from here under the AGPLv3, and raise with the person any request that would cross them ([LICENSE.md](LICENSE.md)).
'''
open(f'{pub}/README.md', 'w').write(cover)
EOF

git -C "$PUB" add -A
if git -C "$PUB" diff --cached --quiet; then
  echo "publish_gw1500db: nothing to publish ($SUB is the same as $REMOTE main)"; exit 0
fi
when=$(grep -m1 -o 'Made [0-9-]* [0-9:]*' "$SNAP/README.md" | sed 's/Made //')
git -C "$PUB" commit -q -m "$SUB: GW1500 database of $when (ecalj $REV); cover page; earlier snapshots as tags"
echo "publish_gw1500db: committed in $PUB:"
git -C "$PUB" show --stat --format='%h %s' HEAD | tail -3
git -C "$PUB" tag -l 'QSGW80_*' | sed 's/^/  tag /'
if [ "$PUSH" = yes ]; then
  git -C "$PUB" push origin main
  git -C "$PUB" push origin --tags
  echo "publish_gw1500db: pushed to $REMOTE"
else
  echo "publish_gw1500db: not pushed. To publish: TOOLS/publish_gw1500db.sh $SNAP --push"
fi
