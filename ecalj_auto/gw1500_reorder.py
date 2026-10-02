#!/usr/bin/env python3
# 2026-10-02 23:25 (user: "以前の結果をサンプルしながら confirm", "水素だけの構造が気になる"): reorder a queue in place.
#   reorder.py <queue> <first mpid,...>   ; run under the queue lock (flock <queue>.lock)
import sys, random, collections
q = sys.argv[1]; first = [m for m in sys.argv[2].split(',') if m]
L = [l for l in open(q).read().split('\n') if l.strip() and not l.startswith('#')]
rec = {l.split('\t')[0]: l for l in L}
cat = lambda m: rec[m].split('\t')[3]; na = lambda m: int(rec[m].split('\t')[2])
head = [m for m in first if m in rec]
rest = [m for m in rec if m not in head]
t1 = [m for m in rest if cat(m) in ('SUSPECT_GOOD', 'MAY_WRONG', 'DRIFT_GOOD', 'UNKNOWN_MAY')]
t2 = [m for m in rest if cat(m) in ('FAILED_MAY', 'NOTCONV_MAY')]
good = [m for m in rest if cat(m) == 'GOOD']
random.seed(20261002)
by = collections.defaultdict(list)
for m in good: by[na(m)].append(m)
for k in by: random.shuffle(by[k])
sample = []                       # round robin over natom: a stratified random order
while any(by.values()):
    for k in sorted(by):
        if by[k]: sample.append(by[k].pop())
out = head + t1
i = j = 0
while i < len(t2) or j < len(sample):
    for _ in range(2):
        if i < len(t2): out.append(t2[i]); i += 1
    if j < len(sample): out.append(sample[j]); j += 1
assert sorted(out) == sorted(rec)
open(q, 'w').write(''.join(rec[m] + '\n' for m in out))
print(len(out), 'lines;', 'head', out[:3])
