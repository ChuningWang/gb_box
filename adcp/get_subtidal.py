import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import netCDF4 as nc
from datetime import datetime

import ttide
from ocean_toolbox import noaa_adcp
import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# stn = 'SEA0847'
stn = 'SEA0850'
# stn = 'SEA1003'
# stn = 'SEA1004'

# -------------- extract data -------------------------------
info = {'stn' : stn,
        'file_dir': out_dir + 'tef/NOAA_ADCP/',
        'sl': 'l',
        'Wp_hrs': 24}

crt = noaa_adcp.get_noaa_current(info)
crt()
time = crt.ctime
U = crt.uraw + 1j*crt.vraw
# u = crt.uraw
# v = crt.vraw
Ubar = U.mean(axis=0)
# ubar = u.mean(axis=0)
# vbar = v.mean(axis=0)
# time convertion
pytime = nc.num2date(time, 'days since 1900-01-01')
yearday = time - \
    nc.date2num(datetime(pytime[0].year, 1, 1), 'days since 1900-01-01') + \
    1
time2 = np.arange(np.floor(time[0]), np.floor(time[0]))
yearday2 = np.arange(np.floor(yearday[0]), np.floor(yearday[-1]))

# -------------- harmonic analysis --------------------------
# perform harmonic analysis
tfit = ttide.t_tide(Ubar, dt=0.1,
                    stime=pytime[0], out_style=None, errcalc='wboot')
m2 = tfit['nameu'] == 'M2  '
ang = tfit['tidecon'][m2, 4][0]-90
Ut = tfit(pytime)

# rotate the speed vector
rot = np.exp(1j*-ang/180*np.pi)
U = U*rot
Ubar = Ubar*rot
Ut = Ut*rot

u = U.real
v = U.imag
ubar = Ubar.real
vbar = Ubar.imag
ut = Ut.real
vt = Ut.imag
utidemax = np.zeros(yearday2.shape)
utidemin = np.zeros(yearday2.shape)
vtidemax = np.zeros(yearday2.shape)
vtidemin = np.zeros(yearday2.shape)
for i, yd in enumerate(yearday2):
    msk = (yearday > yd) & (yearday < yd+1)
    utidemax[i] = ut[msk].max()
    utidemin[i] = ut[msk].min()
    vtidemax[i] = vt[msk].max()
    vtidemin[i] = vt[msk].min()

# -------------- calculate inflow, outflow ------------------
# interpolate
zlevs = 200
z = crt.z
z2 = np.linspace(z[0], z[-1], zlevs)
dz = z2[1]-z2[0]
ures = interp1d(z, crt.u, axis=0)(z2)
vres = interp1d(z, crt.v, axis=0)(z2)

quin = np.zeros(time.shape)
quout = np.zeros(time.shape)
qvin = np.zeros(time.shape)
qvout = np.zeros(time.shape)

for i, t in enumerate(time):
    msk = ures[:, i] >= 0
    quin[i] = ures[msk, i].sum()*dz
    quout[i] = ures[~msk, i].sum()*dz
    msk = vres[:, i] >= 0
    qvin[i] = vres[msk, i].sum()*dz
    qvout[i] = vres[~msk, i].sum()*dz

fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.set_xlim(yearday2[0]+2, yearday2[-1]-1)
ax.set_xlabel(r'Yearday')
ax.set_ylabel(r'Q [m$^2$s$^{-1}$]')
ax2.set_ylabel(r'Tidal Amplitude [m$\cdot$s$^{-1}$]')
# ax.tick_params(axis='y', colors='blue')

ax2.plot(yearday2+0.5, vtidemax, 'k')
# ax.plot(yearday, quin, 'r')
# ax.plot(yearday, quout, '--r')
ax.plot(yearday, qvin, 'b')
ax.plot(yearday, qvout, '--b')

plt.show()

