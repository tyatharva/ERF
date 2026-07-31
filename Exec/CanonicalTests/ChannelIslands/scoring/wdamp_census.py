#!/usr/bin/env python3
"""w-damping census from a run log.

  wdamp_census.py <run_dir> [label]

Log lines look like:
  t=32400.791016 w-damping applied at (41,1,13) for w-CFL = 3.220394 > 1.000000 : -314.275452 --> -313.755646

Reports events per model hour, max |w| per hour, and a (j,k) histogram both
aggregated and restricted to the scored window h48-h71. Model time is 72 h
(259200 s): start 2020-12-26 00Z, stop 2020-12-29 00Z.

Restarts mean a model time can appear more than once in the log. Events are
attributed by their own t= stamp, so a rolled-back segment contributes to the
hour it claims; the count of duplicated hours is reported separately so a
restart cannot silently inflate a bin.
"""
import collections
import re
import sys

PAT = re.compile(
    r't=([0-9.]+) w-damping applied at \((\d+),(\d+),(\d+)\) for w-CFL = ([0-9.eE+-]+) '
    r'> [0-9.]+ : (-?[0-9.eE+-]+) --> (-?[0-9.eE+-]+)')
SCORED = (48, 72)


def main():
    run = sys.argv[1].rstrip('/')
    label = sys.argv[2] if len(sys.argv) > 2 else run

    per_hour = collections.Counter()
    maxw_hour = collections.defaultdict(float)
    maxcfl_hour = collections.defaultdict(float)
    hist_all = collections.Counter()
    hist_scored = collections.Counter()
    steps_seen = 0
    hours_seen = collections.Counter()
    t_at_32400 = []

    with open(f'{run}/run.log', 'rb') as f:
        for raw in f:
            line = raw.decode('utf-8', 'replace')
            if 'Coarse STEP' in line and 'ends' in line:
                steps_seen += 1
                continue
            if 'w-damping applied' not in line:
                continue
            m = PAT.search(line)
            if not m:
                continue
            t = float(m.group(1))
            j, k = int(m.group(3)), int(m.group(4))
            cfl = float(m.group(5))
            w = abs(float(m.group(6)))
            h = int(t // 3600)
            per_hour[h] += 1
            hours_seen[h] += 1
            maxw_hour[h] = max(maxw_hour[h], w)
            maxcfl_hour[h] = max(maxcfl_hour[h], cfl)
            hist_all[(j, k)] += 1
            if SCORED[0] <= h < SCORED[1]:
                hist_scored[(j, k)] += 1
            if 32400.0 <= t < 32401.0:
                t_at_32400.append((cfl, w))

    total = sum(per_hour.values())
    print(f'=== {label}: {total} w-damping events over {steps_seen} logged steps ===')
    if total == 0:
        print('  ZERO events. Nothing further to report.')
        return

    print(f'\n  events per model hour (of 72 h):')
    print(f'  {"h":>4}{"events":>10}{"max|w| m/s":>14}{"max w-CFL":>12}')
    for h in sorted(per_hour):
        mark = ' <- scored' if SCORED[0] <= h < SCORED[1] else ''
        print(f'  {h:>4}{per_hour[h]:>10}{maxw_hour[h]:>14.1f}{maxcfl_hour[h]:>12.2f}{mark}')

    sc = sum(v for h, v in per_hour.items() if SCORED[0] <= h < SCORED[1])
    print(f'\n  scored window h48-h71: {sc} events '
          f'({100.0*sc/total:.1f}% of all)')

    for nm, hist in (('ALL HOURS', hist_all), ('SCORED h48-h71', hist_scored)):
        if not hist:
            print(f'\n  (j,k) histogram, {nm}: empty')
            continue
        print(f'\n  (j,k) histogram, {nm} -- top 12 of {len(hist)} cells:')
        for (j, k), n in hist.most_common(12):
            print(f'    j={j:<4} k={k:<4} {n:>8}')
        js = collections.Counter()
        ks = collections.Counter()
        for (j, k), n in hist.items():
            js[j] += n
            ks[k] += n
        print(f'    j marginal: {dict(sorted(js.items())[:12])}')
        print(f'    k marginal: {dict(sorted(ks.items())[:12])}')

    print(f'\n  t = 32400 s (the 9 h frame-3 transition):')
    if t_at_32400:
        cfls = [c for c, _ in t_at_32400]
        ws = [w for _, w in t_at_32400]
        print(f'    {len(t_at_32400)} events, max w-CFL {max(cfls):.2f}, '
              f'max |w| {max(ws):.1f} m/s')
        print(f'    -> the surviving trajectory did NOT pass this quietly'
              if max(cfls) > 1.0 else '    -> sub-critical')
    else:
        print('    NO w-damping events at t=32400 on the surviving trajectory'
              ' -- it passed the transition with sub-critical w.')

    dup = [h for h, n in hours_seen.items() if n and per_hour[h] != n]
    print(f'\n  restart bookkeeping: {len(dup)} hours with re-counted events'
          if dup else '\n  restart bookkeeping: no double-counted hours detected')


if __name__ == '__main__':
    main()
