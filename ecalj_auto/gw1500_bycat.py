#!/usr/bin/env python3
"""gw1500_bycat.py: GW1500_by_category.md, the 1546 materials listed by category (2026-10-02, user: "GW1500 は分類で見れるように").

    python3 gw1500_bycat.py          (in ecalj_auto/)

Input: gw1500_notes_20261001.tsv (category, notes, May and rerun results, MP data; made by gw1500_notes.py) and
gw1500_status_20260930.tsv (the LDA gap of the May runs). The categories are those of GW1500_status.md table 6.
The adopted gap: the rerun (fp32, 2026-09-30 ~) when it converged (0 for a metal), else the May value when May was GOOD, else none.
Regenerate after either table changes; do not edit the .md by hand.
"""
import csv
from pathlib import Path
D = Path(__file__).parent
notes = list(csv.DictReader(open(D/'gw1500_notes_20261001.tsv'), delimiter='\t'))
lda = {r['mpid']: r['gap_LDA_eV'] for r in csv.DictReader(open(D/'gw1500_status_20260930.tsv'), delimiter='\t')}
ORDER = [
 ('INVALID_STRUCTURE', 'Rb8 だけ（Rb-IV の高圧相を MP が圧力ゼロで緩和、Rb8 の環が離れて並ぶ）。別の化合物の副格子を抜き出した 11 件は 2026-10-06 から 5 月のカテゴリに戻し、注記の「構造:」に書く（user「Rb8 だけ INVALID_STRUCTURE にしておこう」）'),
 ('SUSPECT_STRUCTURE', 'ICSD の備考には出ないが、非整合層状化合物の副格子と同じ形（1 原子の体積が岩塩型の約 2〜4 倍、配位 4〜5）'),
 ('MAY_WRONG', '5 月のギャップが誤り（LDA より小さい）。回し直しの値を使う'),
 ('FAILED_MAY', '5 月に使える結果が無かった。回し直しで全部収束'),
 ('UNKNOWN_MAY', 'LDA でもギャップが無いか 0.1 eV 程度（金属・半金属）'),
 ('SUSPECT_GOOD', '5 月は GOOD だが、QSGW80 < LDA または 2 反復目以降に 0.5 eV を超えて振動。回し直しで全部収束'),
 ('NOTCONV_MAY', '5 月は 10 反復で収束しなかった。回し直しで全部収束'),
 ('DRIFT_GOOD', '5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに 0.08 eV 以上動いた。回し直していない'),
 ('GOOD', '5 月に収束、注記なし。回し直していない'),
]
def adopted(r):
    if r['rerun_verdict'] == 'CONVERGED' and r['rerun_gap_fp32_eV']: return r['rerun_gap_fp32_eV'], '回し直し'
    if r['rerun_verdict'] == 'CONVERGED_METAL': return '0（金属）', '回し直し'
    if r['may_state'] == 'GOOD' and r['may_gap_TF32_eV']: return r['may_gap_TF32_eV'], '5 月'
    return '', ''
link = lambda m: f'[{m}](https://next-gen.materialsproject.org/materials/{m})'
esc = lambda s: s.replace('|', '\\|').replace('\n', ' ')
out = ['# GW1500 の分類ごとの一覧', '',
       f'`gw1500_bycat.py` が [`gw1500_notes_20261001.tsv`](gw1500_notes_20261001.tsv) から作る（手で直さない）。分類の意味と経過は '
       '[GW1500_status.md](GW1500_status.md) §5 の表 6。採用のギャップ: 回し直し（fp32、2026-09-30〜）が収束していればその値、'
       '金属として収束したものは 0、でなければ 5 月が GOOD なら 5 月の値（TF32）。QSGW80（`scaledsigma = 0.8`）、eV。mpid は Materials Project へのリンク。', '',
       '**表 1**. 分類の数', '', '| 分類 | 数 | 意味 |', '| --- | --- | --- |']
groups = {c: [r for r in notes if r['category'] == c] for c, _ in ORDER}
assert sum(len(v) for v in groups.values()) == len(notes), 'a category not in ORDER'
for c, m in ORDER: out.append(f'| [`{c}`](#{c.lower().replace("_", "-")}) | {len(groups[c])} | {m} |')
out.append(f'| 計 | {len(notes)} | |')
for i, (c, m) in enumerate(ORDER, 2):
    rows = sorted(groups[c], key=lambda r: (int(r['natom']), r['formula']))
    out += ['', f'## {c}', '', f'{m}。{len(rows)} 物質。', '']
    if c in ('INVALID_STRUCTURE', 'SUSPECT_STRUCTURE'):
        out += [f'**表 {i}**. 1 原子の体積（Å³）・凸包からの距離（eV/原子）は MP の値。ICSD の備考は MP の項目のもの', '',
                '| mpid | 組成 | 原子 | 1 原子の体積 | 凸包から | 5 月 | 回し直し | ICSD の備考 | 特徴 |', '| --- | --- | --- | --- | --- | --- | --- | --- | --- |']
        for r in rows:
            may = f"{r['may_state']} {r['may_gap_TF32_eV']}".strip(); rr = f"{r['rerun_verdict']} {r['rerun_gap_fp32_eV']}".strip()
            out.append(f"| {link(r['mpid'])} | {r['formula']} | {r['natom']} | {r['mp_vol_per_atom_A3']} | {r['mp_ehull_eV_atom']} | {may} | {rr} | "
                       f"{esc(r['mp_icsd_remarks'])} | {esc(r['note'])} |")
    else:
        out += [f'**表 {i}**', '', '| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |', '| --- | --- | --- | --- | --- | --- | --- |']
        for r in rows:
            g, src = adopted(r)
            out.append(f"| {link(r['mpid'])} | {r['formula']} | {r['natom']} | {g} | {src} | {r['rerun_gapLDA_eV'] or lda.get(r['mpid'], '')} | {esc(r['note'])} |")
(D/'GW1500_by_category.md').write_text('\n'.join(out) + '\n')
print('wrote GW1500_by_category.md:', {c: len(v) for c, v in groups.items()})
