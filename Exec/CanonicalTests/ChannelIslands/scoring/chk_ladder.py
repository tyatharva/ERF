import numpy as np
m=np.load('/app/ERF/mrms_20201228_on_grid.npy')
terr=np.load('/app/ERF/scoring_ab_dav/terrain.npy')
NX,NY=m.shape
ii,jj=np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
d=np.minimum.reduce([ii,jj,NX-1-ii,NY-1-jj])
land=terr>terr.min()+20
masks=[('interior d>=20 (all)',(d>=20)),('interior d>=20 LAND',(d>=20)&land),
       ('LAND all d',land),('whole domain',np.ones_like(land,bool))]
print(f'{"mask":22s} {"n":>6} ' + ' '.join(f'{t:>13}' for t in ('>=1mm','>=5mm','>=15mm','>=30mm')))
for nm,msk in masks:
    a=m[msk&np.isfinite(m)]
    row=f'{nm:22s} {a.size:6d} '
    for thr in (1.,5.,15.,30.):
        n=int((a>=thr).sum()); row+=f'{n:6d} ({n/a.size*100:4.1f}%) '
    print(row)
print('\nFSS needs a workable base rate; a rung with <~100 cells or <2% coverage is not scoreable.')
