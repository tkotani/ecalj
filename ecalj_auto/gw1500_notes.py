#!/usr/bin/env python3
"""GW1500: one row per material with a category and a note on what is wrong or special with it (2026-10-01).

user 2026-10-01: "failed, undecided, suspect, NOTCONV: attach comments; write their features, and the difference of the
conditions". Nothing is removed from the set; the note says what the material is and how its number was obtained.

Inputs (all in this directory):
  gw1500_status_20260930.tsv     the May runs (TF32), GW1500_status.md section 7
  gw1500_mp_20261001.json        Materials Project: hull energy, density, volume, ICSD remarks (read by the API, no key)
  gw1500_rerun_logs_20261001.txt the logs of the reruns on kt1 (run1, try, run2m, run3), '<log> <date> <time> <worker> <mpid> ...'
Output: gw1500_notes_20261001.tsv

  python3 gw1500_notes.py      (again when more rerun logs come; copy the logs of kt1 into gw1500_rerun_logs_20261001.txt)
"""
import csv, json, re

STATUS = 'gw1500_status_20260930.tsv'
MPJSON = 'gw1500_mp_20261001.json'
LOGS = 'gw1500_rerun_logs_20261001.txt'
OUT = 'gw1500_notes_20261001.tsv'

# Structures that are not a material (checked 2026-10-01 by the MP API, the POSCARs and the ICSD remarks).
# Most are one sublattice cut out of another compound (misfit layer compounds, hydroxides, intercalated graphite):
# isolated, it gets a molecular HOMO-LUMO gap, which is how it passed the selection (PBE gap > 0).
# 2026-10-06 14:56, user: "Rb8 だけ INVALID_STRUCTURE にしておこう": only Rb8 keeps the category; the sublattices below keep their
# note (SUBLATTICE) and the category of their May runs.
INVALID = {
 'mp-1179832': 'Rb-IV（高圧相、ICSD 109016）を MP が圧力ゼロで緩和。一辺 19.9 Å の胞に Rb8 の正八角形の環（Rb–Rb 4.63 Å、隣 2 個）。1 原子 551 Å³（bcc Rb 93）、凸包から 0.45 eV/原子。ギャップは環の HOMO–LUMO。平面波 3.2 万で lmf --jobgw=1 が 1 ランク 35 GB',
}
SUBLATTICE = {
 'mp-1056418': 'Sr–Co–O の「(Sr) part」を抜き出したもの。10.99×10.99×4.42 Å の胞に Sr 1 個（Sr の鎖、隣 2 個、4.42 Å）。1 原子 462 Å³（fcc Sr 56）、凸包から 1.40 eV/原子',
 'mp-730101': 'NH4D2PO4 の H（D）だけを抜き出したもの。H2 分子 4 個（H–H 0.74 Å）、密度 0.04 g/cm³',
 'mp-554134': 'Gd–Sn–Nb–S の非整合層状化合物の SnS 層。a=4.13、c=22.4 Å、1 原子 95.7 Å³（SnS 約 24）、配位 4',
 'mp-726184': 'Pb–Ti–S の非整合層状化合物の PbS 層。a=4.22、c=16.6 Å、1 原子 73.6 Å³（岩塩型 PbS 26.7、最近接 2.99 Å・配位 6 に対し 2.66 Å・配位 5）',
 'mp-727323': 'Pb–Ti–S の非整合層状化合物の「Pb S-part」。a=4.22、c=11.2 Å、1 原子 49.7 Å³（岩塩型 PbS 26.7）',
 'mp-727322': 'Sn–Ti–S の非整合層状化合物の SnS 層。a=4.09、c=11.0 Å、1 原子 45.3 Å³',
 'mp-8781': 'Sn–Nb–S の非整合層状化合物の SnS 層。a=4.12、c=11.0 Å、1 原子 46.8 Å³',
 'mp-1062030': 'Pb–Ti–S の非整合層状化合物の「Ti S2-part」。c=11.25 Å（1T-TiS2 は 5.70）、1 原子 37.9 Å³（約 19）。5 月は 1 反復目に 0.65 → 2.49 eV と跳んだ',
 'mp-1079707': 'Ca–Co 水酸化物の「(Ca(OH))-part」から H を除いたもの（Ca4O4）。c=16.75 Å の層、密度 2.00（岩塩型 CaO 3.34）。5 月は 9.97 → 3.31 → 4.78 eV と荒れた',
 'mp-569304': '硝酸を挿入した黒鉛（Graphite, nitrated）から N・O を除いたもの。密度 1.38（黒鉛 2.26）、層間が開いたまま',
 'mp-569416': '硝酸を挿入した黒鉛から N・O を除いたもの。密度 1.67（黒鉛 2.26）',
}
# 2026-10-07 14:54: the Sr chain cut out of Sr-Co-O; its MLO model is skipped, the QSGW80 value kept. Category TOO_THEORETICAL until
# 2026-10-07 18:51, then TOO_LARGE like Rb8 (user: "どっちも TOO_Large かな"; 462 A^3 per atom)
TOO_THEORETICAL = {'mp-1056418'}
# The same shape as the misfit sublattices above (a ~ 4.2 A, c 11-22 A, rock-salt-like layers with coordination 4-5),
# though the ICSD remarks do not name the other compound.
SUSPECT_STRUCT = {
 'mp-561320': 'PbS。a=4.22、c=22.3 Å、1 原子 98.9 Å³（岩塩型 26.7）、配位 5。非整合層状化合物の PbS 層と同じ形',
 'mp-20526': 'PbS。a=4.22、c=11.3 Å、1 原子 50.2 Å³、配位 4。非整合層状化合物の PbS 層と同じ形',
 'mp-22009': 'PbSe。1 原子 58.5 Å³、配位 5（岩塩型 PbSe は約 29.5、配位 6）。層状化合物の副格子と見られる',
 'mp-8936': 'SnSe。a=4.29、c=11.0 Å、1 原子 50.5 Å³、配位 5。非整合層状化合物の SnSe 層と同じ形',
}
MAY_WRONG = {'mp-546711': '5 月は GOOD でギャップ 0.26 eV（LDA 5.41 より小さい）。回し直しで 8.95 eV。5 月の値は誤り'}

COND_RERUN = 'fp32; 989a18637 (2026-09-29); POSCAR→vasp2ctrl→ctrlgenToml --ssig=0.8 (2026-09 の雛形); t_tetrakbt=300 t_sigmaw=300; LDA から; 最大 10 反復 (NaCoO2 20、金属 15)'


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def may_conditions(r):
    s = f"TF32; 本計算 {r['prod_input']} {r['prod_bindir']}"
    n = int(r['n_ac_runs'] or 0)
    if n:
        s += f" → 追加の計算 TOML(05-11) ×{n}"
    if r['redo_status']:
        s += ' → 最初から REDO(TOML)'
    return s + '; χ0 T=0, Σ GaussSmear esmr=0.003 Ry'


def read_logs():
    res = {}
    rx = re.compile(r'^(\S+) (\S+ \S+) (\S+) (mp-\d+) (\S+) iter=(\d+) gapLDA=(\S+) gap=(\S+) (\d+)s(?: dqp=(\S+))?')
    for line in open(LOGS):
        m = rx.match(line)
        if not m:
            continue
        log, dt, wk, mp, v, it, gl, g, sec, dqp = m.groups()
        place = log.split('/')[0]
        if 'continued from iter' in line:
            place += '(続き)'
        res[mp] = dict(place=place, verdict=v, iter=it, gapLDA=gl, gap=g, sec=sec, dqp=dqp or '', when=dt)  # the last line wins
    return res


def main():
    st = list(csv.DictReader(open(STATUS), delimiter='\t'))
    mp = json.load(open(MPJSON))
    rr = read_logs()
    q = {}
    try:  # run3 priorities (the reason a GOOD of May was rerun)
        for line in open('gw1500_recheck_queue_20260930.tsv'):
            if not line.startswith('#'):
                w = line.rstrip('\n').split('\t')
                q[w[0]] = w[5]
    except OSError:
        pass
    rows = []
    for r in st:
        m = r['mpid']; x = mp.get(m, {}); R = rr.get(m)
        notes = []
        if m in TOO_THEORETICAL:   # 2026-10-07 14:54 (user: "Sr は skipped (too theoretical) とかにする？")
            cat = 'TOO_LARGE'; notes.append('胞が大きい（1 原子 462 Å³）。QSGW80 の値は残し、MLO は作らない（skipped。2026-10-07 18:51 まで TOO_THEORETICAL）')   # the structure: SUBLATTICE below
        elif m in INVALID:
            cat = 'TOO_LARGE'; notes.append(INVALID[m])   # INVALID_STRUCTURE until 2026-10-06 19:28 (user: "大きすぎる、というべきかな")
        elif m in SUSPECT_STRUCT:
            cat = 'SUSPECT_STRUCTURE'; notes.append(SUSPECT_STRUCT[m])
        elif m in MAY_WRONG:
            cat = 'MAY_WRONG'; notes.append(MAY_WRONG[m])
        elif r['final_state'] == 'FAILED':
            cat = 'FAILED_MAY'; notes.append(f"5 月の失敗: {r['final_detail']}")
        elif r['final_state'] == 'UNKNOWN':
            cat = 'UNKNOWN_MAY'; notes.append(f"5 月: {r['final_detail']}（金属・半金属。gwscconv がギャップを読めず止まった）")
        elif q.get(m, '').startswith('1:'):
            cat = 'SUSPECT_GOOD'; why = q[m][2:].replace('SUSPECT:', '')
            notes.append(f"5 月は GOOD だが {'QSGW80 のギャップが LDA より小さい' if 'gap<LDA' in why else '2 反復目以降に 0.5 eV を超えて振動'}（履歴 {r['gap_history_eV']}）")
        elif r['final_state'] == 'NOTCONV':
            cat = 'NOTCONV_MAY'; notes.append(f"5 月は 10 反復で収束せず（最後の 3 反復の幅 {r['last3_swing_eV']} eV、履歴 {r['gap_history_eV']}）")
        elif q.get(m, '').startswith('3:'):
            cat = 'DRIFT_GOOD'; notes.append(f"5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに {q[m].split(':')[-1]} eV 動いていた（履歴 {r['gap_history_eV']}）。回し直していない")
        else:
            cat = 'GOOD'
        if not x.get('formula_pretty'):
            notes.append('MP の summary が無い（ID が消えたか統合された）')
        if m in SUBLATTICE:
            notes.insert(0, '構造: ' + SUBLATTICE[m])
        if R:
            g, gm = f(R['gap']), f(r['gap_QSGW80_last_llmf_eV'])
            s = f"回し直し({R['place']}): {R['verdict']} {R['iter']} 反復、ギャップ {R['gap']} eV（LDA {R['gapLDA']}）"
            if g is not None and gm is not None:
                s += f"、5 月との差 {g - gm:+.2f} eV"
            if cat == 'TOO_LARGE':
                s += '。構造が不正なので値に物理的な意味は無い'
            notes.append(s)
        elif cat in ('NOTCONV_MAY',):
            notes.append('回し直し待ち（run3）')
        va = (x['volume'] / x['nsites']) if x.get('volume') and x.get('nsites') else None
        rem = '; '.join(t for t in (x.get('remarks') or []) if t)[:120]
        rows.append({
            'mpid': m, 'formula': r['formula'], 'natom': r['natom'], 'category': cat, 'note': ' / '.join(notes),
            'may_state': r['final_state'], 'may_gap_TF32_eV': r['gap_QSGW80_last_llmf_eV'], 'may_conditions': may_conditions(r),
            'rerun_place': R['place'] if R else '', 'rerun_verdict': R['verdict'] if R else '', 'rerun_iter': R['iter'] if R else '',
            'rerun_gap_fp32_eV': R['gap'] if R else '', 'rerun_gapLDA_eV': R['gapLDA'] if R else '',
            'rerun_conditions': COND_RERUN if R else '',
            'mp_ehull_eV_atom': f"{x['energy_above_hull']:.3f}" if x.get('energy_above_hull') is not None else '',
            'mp_density_gcm3': f"{x['density']:.2f}" if x.get('density') else '', 'mp_vol_per_atom_A3': f"{va:.1f}" if va else '',
            'mp_icsd_remarks': rem})
    with open(OUT, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    from collections import Counter
    c = Counter(r['category'] for r in rows)
    print(OUT, len(rows), dict(c))
    for cat in ('FAILED_MAY', 'UNKNOWN_MAY', 'SUSPECT_GOOD', 'NOTCONV_MAY', 'TOO_LARGE', 'TOO_THEORETICAL', 'SUSPECT_STRUCTURE'):
        sub = [r for r in rows if r['category'] == cat]
        print(f"{cat:18s} {len(sub):4d}  rerun: {dict(Counter(r['rerun_verdict'] or 'not yet' for r in sub))}")


if __name__ == '__main__':
    main()
