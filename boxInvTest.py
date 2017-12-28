""" test gb_box. """

from datetime import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from gb_box import boxPrep
from gb_box import boxODE

INFO = {'case': 'Box2',
        'box_method': 2,
        'clim': True,
        'compare': True,
        't0': nc.date2num(datetime(2008, 01, 01), 'days since 1900-01-01'),
        't1': nc.date2num(datetime(2009, 01, 01), 'days since 1900-01-01'),
        'grid_file_name': '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_lr_grd.nc',
        'river_file_name': '/Users/CnWang/git/gb_box/data/gb_box_rivers2.nc',
        'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
        'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc'}

box = boxPrep.Box(INFO)
box()

ms = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366])

f0data = box.BoxRivers.climatology['river0']*box.consts['s2d']
ftime = box.BoxRivers.climatology['time']
f0 = np.zeros(12)
for i in range(12):
    msk = (ftime >= ms[i]) & (ftime < ms[i+1])
    f0[i] = f0data[msk].mean()

time = box.measurements['time']
Vu0 = box.volumes['Vu0']
Vi0 = box.volumes['Vi0']
dSu0 = box.measurements['dSu0']
dSi0 = box.measurements['dSi0']
Su0 = box.measurements['Su0']
Si0 = box.measurements['Si0']
Su2 = box.measurements['Su2']
Si2 = box.measurements['Si2']

q0 = np.zeros(12)
w0 = np.zeros(12)
for i in range(12):
    B = np.array([Vu0*dSu0[i]+f0[i]*Si0[i],
                  Vi0*dSi0[i]+f0[i]*(Si2[i]-Si0[i])])
    A = np.array([[Si0[i]-Su0[i], Si0[i]-Su0[i]],
                  [Si2[i]-Si0[i], Su0[i]-Si0[i]]])
    slv = np.linalg.solve(A, B)
    q0[i] = slv[0]
    w0[i] = slv[1]


# q0 = (Vu0*dSu0+Si0*f0)/(Si0-Su0)

fig, (ax0, ax1, ax2) = plt.subplots(3, sharex=True)
ax0.set_xlim(0, 366)
ax0.plot(time, Su0)
ax0.plot(time, Si0)
ax0.plot(time, Si2, 'r')
ax1.plot(time, f0)
ax1.plot(time, q0)
ax2.plot(time, w0)

fig.savefig('inv.png', dpi=500)
plt.close()
