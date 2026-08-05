import numpy as np, netCDF4 as nc
f=nc.Dataset('/app/ERF/wrfout_d02_2020-12-28_00_00_00')
for a in ('BUCKET_MM','BUCKET_J','PREC_ACC_DT'):
    if hasattr(f,a): print(f'  global {a} = {getattr(f,a)}')
R0=np.asarray(f.variables['RAINNC'][0]);  R23=np.asarray(f.variables['RAINNC'][-1])
B0=np.asarray(f.variables['I_RAINNC'][0]);B23=np.asarray(f.variables['I_RAINNC'][-1])
B=float(getattr(f,'BUCKET_MM',100.0))
naive=R23-R0
tot=(R23+B*B23)-(R0+B*B0)
wrapped=int((B23>B0).sum())
print(f'\nbucket_mm used = {B}')
print(f'cells whose bucket incremented in the window: {wrapped} of {R0.size} ({wrapped/R0.size*100:.2f}%)')
print(f'naive RAINNC diff : min {naive.min():9.3f}  max {naive.max():8.2f}  mean {naive.mean():7.3f}  n_negative {int((naive<0).sum())}')
print(f'bucket-corrected  : min {tot.min():9.3f}  max {tot.max():8.2f}  mean {tot.mean():7.3f}  n_negative {int((tot<0).sum())}')
print(f'RAINC over window : {float((np.asarray(f.variables["RAINC"][-1])-np.asarray(f.variables["RAINC"][0])).max()):.4f} mm max (0 => no cumulus scheme, expected at 1.5 km)')
np.save('/app/ERF/wrf_d02_accum_mm.npy', tot)
np.save('/app/ERF/wrf_d02_latlon.npy',
        np.stack([np.asarray(f.variables['XLAT'][0]), np.asarray(f.variables['XLONG'][0])]))
print('saved accumulation + coords')
