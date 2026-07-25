import yt, numpy as np, glob
yt.set_log_level(50)
f={}
for cfg in ['p1_base','p1_rho','p1_rhow']:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    f[cfg]=(np.asarray(g[('boxlib','density')]), np.asarray(g[('boxlib','z_velocity')]))
for c in ['p1_rho','p1_rhow']:
    dr=f[c][0]-f['p1_base'][0]; dw=f[c][1]-f['p1_base'][1]
    print(f'{c} vs base: max|drho| {np.abs(dr).max():.6e}  max|dw| {np.abs(dw).max():.6e}  '
          f'outermost-row |drho| {np.abs(dr[0,:,:]).max():.6e}')
