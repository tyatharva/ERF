#!/usr/bin/env python3
"""Test 2: MODE-style object verification -- the standard answer to neighbourhood
scores rewarding smoothness.

Objects are defined the MODE way: convolve, threshold, label connected
components, discard specks. Forecast objects are matched to observed objects by
an interest score combining centroid separation and area similarity. Reported:
matched fraction, median centroid displacement, area bias, intensity bias.

Unlike FSS this never compares fields cell-by-cell, so a smooth field earns
nothing from being smooth -- it either puts an object near an observed one or it
does not.

Control: a field scored against ITSELF must give 100% matched, 0 km displacement,
and area/intensity bias 1.0.
"""
import sys
import numpy as np
from scipy.ndimage import gaussian_filter, label, center_of_mass

NX, NY, DX = 192, 96, 3.0
CONV_SIGMA = 1.0        # MODE convolution radius, cells
MIN_AREA = 10           # cells; discard specks
D0 = 100.0              # km, interest half-scale


def objects(field, thr, mask):
    sm = gaussian_filter(np.nan_to_num(field), CONV_SIGMA)
    binm = (sm >= thr) & mask
    lab, n = label(binm, structure=np.ones((3, 3)))
    out = []
    for i in range(1, n + 1):
        sel = lab == i
        a = int(sel.sum())
        if a < MIN_AREA:
            continue
        cy, cx = center_of_mass(sel)
        vals = np.nan_to_num(field)[sel]
        out.append(dict(area=a, cx=cy * DX, cy=cx * DX,
                        mean=float(vals.mean()), p90=float(np.percentile(vals, 90))))
    return out


def match(fo, oo):
    pairs = []
    used = set()
    for o in oo:
        best, bi = -1.0, None
        for k, f in enumerate(fo):
            if k in used:
                continue
            dist = np.hypot(f['cx'] - o['cx'], f['cy'] - o['cy'])
            ar = min(f['area'], o['area']) / max(f['area'], o['area'])
            interest = 0.5 * max(0.0, 1.0 - dist / D0) + 0.5 * ar
            if interest > best:
                best, bi = interest, k
        if bi is not None and best > 0.25:
            used.add(bi)
            pairs.append((fo[bi], o, np.hypot(fo[bi]['cx'] - o['cx'], fo[bi]['cy'] - o['cy'])))
    return pairs


def main():
    arms = [a.split('=', 1) for a in sys.argv[1:]]
    base = arms[0][1].rstrip('/') + '/'
    mrms = np.load(base + 'mrms_mm.npy')
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    mask = (d >= 20) & np.isfinite(mrms)

    fields = {}
    for lab_, dd in arms:
        dd = dd.rstrip('/') + '/'
        fields[lab_] = np.nan_to_num(
            np.load(dd + ('era5_mm.npy' if lab_.startswith('ERA5') else 'erf_mm.npy')))

    for thr in (5.0, 15.0, 30.0):
        oo = objects(mrms, thr, mask)
        ctl = match(objects(mrms, thr, mask), oo)
        assert len(ctl) == len(oo) and max(p[2] for p in ctl) < 1e-9, 'self-match control failed'
        print(f'\n=== objects at {thr:.0f} mm (conv sigma {CONV_SIGMA} cells, '
              f'min area {MIN_AREA} cells), interior d>=20 ===')
        print(f'  MRMS: {len(oo)} objects, total area {sum(o["area"] for o in oo)} cells; '
              f'self-match control OK (0 km, 100%)')
        print(f'{"arm":<14} {"n_obj":>6} {"matched":>8} {"med displ km":>13} '
              f'{"area bias":>10} {"intens bias":>12}')
        for lab_, f in fields.items():
            fo = objects(f, thr, mask)
            pairs = match(fo, oo)
            if not pairs:
                print(f'{lab_:<14} {len(fo):>6} {"0/%d" % len(oo):>8}   no matches')
                continue
            disp = np.median([p[2] for p in pairs])
            ab = np.median([p[0]['area'] / p[1]['area'] for p in pairs])
            ib = np.median([p[0]['p90'] / p[1]['p90'] for p in pairs if p[1]['p90'] > 0])
            print(f'{lab_:<14} {len(fo):>6} {f"{len(pairs)}/{len(oo)}":>8} '
                  f'{disp:>13.1f} {ab:>10.2f} {ib:>12.2f}')


if __name__ == '__main__':
    main()
