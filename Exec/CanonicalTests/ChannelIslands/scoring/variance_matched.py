#!/usr/bin/env python3
"""Test 1 + 3: is Davies' sub-200 km placement advantage skill, or smoothness?

Percentile matching removes AMPLITUDE bias but not VARIANCE bias. Davies carries
1/150th of NSCBC's 8-19 km spectral power, and binarising a near-smooth field at
its own quantile yields large contiguous blobs that overlap a structured
observation field more readily than a structured forecast does.

Method: smooth each structured field with a Gaussian until its 8-19 km spectral
ratio against MRMS MATCHES Davies' 0.003, then re-score percentile-matched FSS.
If the advantage collapses, the sub-200 km gap was variance, not placement.

ERA5 is carried through as the control -- it is also smooth relative to MRMS, so
if its advantage moves the same way under the same operation, the metric is
implicated rather than any particular model.

Controls: sigma=0 must reproduce the unsmoothed scores exactly; the spectral
ratio must be monotone in sigma so the match is well defined.
"""
import sys
import numpy as np
from scipy.ndimage import uniform_filter, gaussian_filter

NX, NY = 192, 96


def spec_ratio(a, ref, mask):
    """8-19 km radial spectral power ratio, same definition as score_metrics."""
    def spec(x):
        x = np.where(mask, np.nan_to_num(x), 0.0)
        core = x[20:NX - 20, 20:NY - 20]
        core = core - core.mean()
        w = np.hanning(core.shape[0])[:, None] * np.hanning(core.shape[1])[None, :]
        F = np.abs(np.fft.rfft2(core * w)) ** 2
        kx = np.fft.fftfreq(core.shape[0], 3.0)[:, None]
        ky = np.fft.rfftfreq(core.shape[1], 3.0)[None, :]
        k = np.sqrt(kx ** 2 + ky ** 2)
        bins = np.linspace(0, k.max(), 40)
        idx = np.digitize(k.ravel(), bins)
        P = np.array([F.ravel()[idx == i].mean() if (idx == i).any() else np.nan
                      for i in range(1, len(bins))])
        lam = 1.0 / (0.5 * (bins[1:] + bins[:-1]))
        return lam, P
    lam, Pa = spec(a)
    _, Pr = spec(ref)
    band = (lam >= 8) & (lam <= 19)
    return float(np.nanmean(Pa[band] / Pr[band]))


def fss(fb, ob, w, mask):
    Pf = uniform_filter(fb.astype(float), size=w, mode='constant')
    Po = uniform_filter(ob.astype(float), size=w, mode='constant')
    num = np.nanmean((Pf[mask] - Po[mask]) ** 2)
    den = np.nanmean(Pf[mask] ** 2) + np.nanmean(Po[mask] ** 2)
    return 1.0 - num / den if den > 0 else np.nan


def binarize(field, mask, rate):
    return field >= float(np.nanquantile(field[mask], 1.0 - rate))


def main():
    arms = [a.split('=', 1) for a in sys.argv[1:]]
    base = arms[0][1].rstrip('/') + '/'
    mrms = np.load(base + 'mrms_mm.npy')
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    mask = (d >= 20) & np.isfinite(mrms)
    valid = np.isfinite(mrms)

    fields = {}
    for lab, dd in arms:
        dd = dd.rstrip('/') + '/'
        fields[lab] = np.nan_to_num(
            np.load(dd + ('era5_mm.npy' if lab.startswith('ERA5') else 'erf_mm.npy')))

    target = spec_ratio(fields['Davies'], mrms, valid)
    print(f'target spectral ratio (Davies) = {target:.4f}')
    for lab in fields:
        print(f'  {lab:<14} native ratio {spec_ratio(fields[lab], mrms, valid):.4f}')

    # find the Gaussian sigma that brings each field to Davies' spectral ratio
    sig = {}
    for lab, f in fields.items():
        r0 = spec_ratio(f, mrms, valid)
        if r0 <= target:
            sig[lab] = 0.0
            continue
        lo, hi = 0.0, 12.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if spec_ratio(gaussian_filter(f, mid), mrms, valid) > target:
                lo = mid
            else:
                hi = mid
        sig[lab] = 0.5 * (lo + hi)
    print('\nGaussian sigma (cells) needed to reach Davies\' variance:')
    for lab in fields:
        sm = gaussian_filter(fields[lab], sig[lab]) if sig[lab] > 0 else fields[lab]
        print(f'  {lab:<14} sigma={sig[lab]:5.2f}  ->  ratio {spec_ratio(sm, mrms, valid):.4f}')

    for phys in (1.0, 5.0):
        rate = float((mrms[mask] >= phys).mean())
        ob = binarize(mrms, mask, rate)
        assert abs(fss(ob, ob, 3, mask) - 1.0) < 1e-12, 'self-FSS control failed'
        print(f'\n=== PM-FSS 3 km, base rate {rate:.4f} (MRMS >= {phys:.0f} mm) ===')
        print(f'{"arm":<14} {"native":>8} {"variance-matched":>18} {"change":>8}')
        for lab, f in fields.items():
            a = fss(binarize(f, mask, rate), ob, 1, mask)
            sm = gaussian_filter(f, sig[lab]) if sig[lab] > 0 else f
            b = fss(binarize(sm, mask, rate), ob, 1, mask)
            if sig[lab] == 0.0:
                assert abs(a - b) < 1e-12, 'sigma=0 identity control failed'
            print(f'{lab:<14} {a:>8.3f} {b:>18.3f} {b-a:>+8.3f}')


if __name__ == '__main__':
    main()
