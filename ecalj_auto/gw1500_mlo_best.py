#!/usr/bin/env python3
"""The variant of the MLO model taken for a material (2026-10-05): best grade PASS > OK > FAIL (grade.json of gw1500_mlo_grade.py);
on a tie the simpler variant (b1, then b2, then b2all), but among FAIL the one with the smallest largest deviation.
    gw1500_mlo_best.py <out of gw1500_mlo_std.sh>     prints "<variant> <grade>"
"""
import json, os, sys
d = sys.argv[1]; R = {'PASS': 0, 'OK': 1, 'FAIL': 2}
c = []
for i, v in enumerate(('b1', 'b2', 'b2all')):
    try:
        g = json.load(open(os.path.join(d, v, 'grade.json')))
    except Exception:
        continue
    c.append((R.get(g['grade'], 3), g['max'] if g['grade'] == 'FAIL' else 0.0, i, v, g['grade']))
print(f'{min(c)[3]} {min(c)[4]}' if c else 'none NONE')
