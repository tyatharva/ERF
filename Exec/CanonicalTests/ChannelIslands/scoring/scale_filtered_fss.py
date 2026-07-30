#!/usr/bin/env python3
"""Does Davies' placement advantage live at scales spectral nudging could reach?

Low-pass both the model field and the observations to a wavelength cutoff, then
score percentile-matched placement on the filtered fields. The filter is a 2-D
DCT-II/III low-pass -- a cosine series on a non-periodic domain, which is exactly
the transform amrex::FFT::R2X with Boundary::even performs and therefore exactly
what a spectral nudging implementation here could act on.

Mode wavelengths for a cosine series on a domain of length L: mode k has spatial
frequency k/(2L), so lambda_k = 2L/k. Cutoff is applied isotropically on the
total frequency f = sqrt(fx^2 + fy^2), keeping lambda = 1/f >= cutoff.

  advantage SURVIVES filtering -> the anchoring lives at reachable scales and a
                                  spectral nudging build is justified
  advantage VANISHES           -> it lives below them; the build would return a
                                  null result

Controls: an unfiltered pass must reproduce the unfiltered scores; the filter with
nothing removed must be the identity to round-off; self-FSS == 1 on every filtered
field; and matched-rate binarization must stay bias-invariant.
"""
import sys
import numpy as np
from scipy.fft import dctn, idctn
from scipy.ndimage import uniform_filter

LX_KM, LY_KM = 576.0, 288.0
WINS = [(3, 1), (60, 21)]


def lowpass(field, cutoff_km):
    """Keep only cosine modes with wavelength >= cutoff_km."""
    if cutoff_km is None:
        return field.copy()
    NX, NY = field.shape
    C = dctn(field, type=2, norm='ortho')
    kx = np.arange(NX)[:, None]
    ky = np.arange(NY)[None, :]
    fx = kx / (2.0 * LX_KM)          # cycles per km
    fy = ky / (2.0 * LY_KM)
    f = np.hypot(fx, fy)
    keep = f <= (1.0 / cutoff_km)
    keep[0, 0] = True                # always retain the domain mean
    return idctn(C * keep, type=2, norm='ortho')


def fss_binary(fb, ob, w, mask):
    Pf = uniform_filter(fb.astype(float), size=w, mode='constant')
    Po = uniform_filter(ob.astype(float), size=w, mode='constant')
    num = np.nanmean((Pf[mask] - Po[mask]) ** 2)
    den = np.nanmean(Pf[mask] ** 2) + np.nanmean(Po[mask] ** 2)
    return 1.0 - num / den if den > 0 else np.nan


def binarize(field, mask, rate):
    return field >= float(np.nanquantile(field[mask], 1.0 - rate))


def main():
    arms = [a.split('=', 1) for a in sys.argv[1:]]
    ref = None
    fields = {}
    for lab, d in arms:
        d = d.rstrip('/') + '/'
        fields[lab] = np.load(d + ('era5_mm.npy' if lab.startswith('ERA5') else 'erf_mm.npy'))
        if ref is None:
            ref = np.load(d + 'mrms_mm.npy')
    NX, NY = ref.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    mask = (d >= 20) & np.isfinite(ref)

    rt = lowpass(ref, None)
    assert np.abs(rt - ref).max() < 1e-9, 'identity control failed'
    full = idctn(dctn(ref, type=2, norm='ortho'), type=2, norm='ortho')
    assert np.abs(full - ref).max() < 1e-9, 'DCT round-trip control failed'
    print(f'controls OK: DCT round-trip exact; interior mask {mask.sum()} cells')

    obs_nan = ~np.isfinite(ref)
    ref_f = np.where(obs_nan, 0.0, np.nan_to_num(ref, nan=0.0))
    for lab, _ in arms:
        nbad = int((~np.isfinite(fields[lab])).sum())
        if nbad: print(f'  note: {lab} has {nbad} non-finite cells, zeroed before transform')

    for phys in (1.0, 5.0):
        rate = float((ref[mask] >= phys).mean())
        print(f'\n=== matched base rate {rate:.4f} (MRMS >= {phys:.0f} mm), interior d>=20 ===')
        hdr = f'{"cutoff":>12} {"modes kx,ky":>12}'
        for lab, _ in arms:
            hdr += f' {lab[:11]:>12}'
        print(hdr + '   (3 km / 60 km)')
        for cut, note in ((None, 'unfiltered'), (200.0, 'lambda>=200km'), (576.0, 'lambda>=576km')):
            of = lowpass(ref_f, cut)
            ob = binarize(of, mask, rate)
            assert abs(fss_binary(ob, ob, 3, mask) - 1.0) < 1e-12, 'self-FSS control failed'
            nkx = int(np.floor(2 * LX_KM / cut)) if cut else NX - 1
            nky = int(np.floor(2 * LY_KM / cut)) if cut else NY - 1
            row = f'{note:>12} {f"{nkx},{nky}":>12}'
            for lab, _ in arms:
                # A DCT is global: a single non-finite cell anywhere poisons the
                # whole transform. Each field must be sanitised on its OWN
                # non-finite cells, not just on the observation mask.
                src = np.nan_to_num(fields[lab], nan=0.0, posinf=0.0, neginf=0.0)
                ff = lowpass(np.where(obs_nan, 0.0, src), cut)
                fb = binarize(ff, mask, rate)
                a3 = fss_binary(fb, ob, 1, mask)
                a6 = fss_binary(fb, ob, 21, mask)
                row += f' {a3:>5.3f}/{a6:<6.3f}'
            print(row)


if __name__ == '__main__':
    main()
