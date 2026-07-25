"""Exact geometry numbers for the trade curve, from geometry_lcc.npz."""
import numpy as np
d = np.load('/app/ERF/precip_check/geometry_lcc.npz')
isl = d['islands']; mx, my = d['mainland_x'], d['mainland_y']
names = ['SantaCruz','SantaRosa','Catalina','SanClemente','SanNicolas',
         'SanMiguel','LosCoronados','SantaBarbaraI','?_sd2','?_la','Anacapa',
         '?_e1','?_sd3','?_s1']
# distance from each island bbox to nearest mainland pixel
print(f"{'island':>14} {'bbox x(km)':>16} {'bbox y(km)':>16} {'dist->mainland(km)':>19}")
for k, (s, xlo, xhi, ylo, yhi, h) in enumerate(isl):
    dx = np.maximum(np.maximum(xlo - mx, mx - xhi), 0.0)
    dy = np.maximum(np.maximum(ylo - my, my - yhi), 0.0)
    dist = np.sqrt(dx**2 + dy**2).min()
    print(f"{names[k]:>14} {xlo/1e3:7.1f}..{xhi/1e3:6.1f} {ylo/1e3:7.1f}..{yhi/1e3:6.1f} {dist/1e3:12.1f}")

# mainland min-x profile vs y (where is the coast?) in 20-km y bands
print("\ncoast profile: min mainland x in y-bands (km):")
for y0 in range(-140, 121, 20):
    m = (my >= y0*1e3) & (my < (y0+20)*1e3)
    if m.sum():
        print(f"  y {y0:4d}..{y0+20:4d}: min x {mx[m].min()/1e3:7.1f}")
# mainland min-y profile vs x (northern coast) in 20-km x bands
print("\ncoast profile: min mainland y in x-bands (km):")
for x0 in range(-200, 201, 25):
    m = (mx >= x0*1e3) & (mx < (x0+25)*1e3)
    if m.sum():
        print(f"  x {x0:4d}..{x0+25:4d}: min y {my[m].min()/1e3:7.1f}")
