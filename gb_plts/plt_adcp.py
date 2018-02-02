from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib import gridspec
import matplotlib as mpl
import netCDF4 as nc
from scipy.signal import filtfilt

from cmocean import cm
from ocean_toolbox import noaa_adcp
import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# -------------- functionals --------------------------------
def filter(data, dt, dt_filter):
    """ filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered

# -------------- read in station data -----------------------
# info = {'stn' : 'SEA1004',
# info = {'stn' : 'SEA0847',
info = {'stn' : 'SEA1010',
        'file_dir': out_dir + 'tef/',
        'sl': 's',
        'Wp_hrs': 2}

crt = noaa_adcp.get_noaa_current(info)
crt()

z = crt.z
time = crt.ctime
pytime = nc.num2date(time, 'days since 1900-01-01')
yr = pytime[0].year
yearday = time - nc.date2num(datetime(yr, 1, 1), 'days since 1900-01-01') + 1
yearday2 = np.arange(np.floor(yearday[0]), np.floor(yearday[-1]))

# ang = 0.35
# ang = np.pi/4
# ang = 0.3
ang = 0
U = crt.u + 1j*crt.v
U = U*np.exp(-ang*1j)
u = U.real
v = U.imag

utbar = u.mean(axis=1)
vtbar = v.mean(axis=1)

u2 = u - np.tile(utbar, (len(yearday), 1)).T
v2 = v - np.tile(vtbar, (len(yearday), 1)).T

utide = u2.mean(axis=0)
vtide = v2.mean(axis=0)

ures = u2 - np.tile(utide, (len(z), 1))
vres = v2 - np.tile(vtide, (len(z), 1))

utidemax = np.zeros(yearday2.shape)
utidemin = np.zeros(yearday2.shape)
vtidemax = np.zeros(yearday2.shape)
vtidemin = np.zeros(yearday2.shape)
for i, yd in enumerate(yearday2):
    msk = (yearday > yd) & (yearday < yd+1)
    utidemax[i] = utide[msk].max()
    utidemin[i] = utide[msk].min()
    vtidemax[i] = vtide[msk].max()
    vtidemin[i] = vtide[msk].min()

# -------------- make plots ---------------------------------
# cmint = -3.
# cmaxt = 3.
# cmine = -1.
# cmaxe = 1.
# cminr = -0.5
# cmaxr = 0.5
# dmin = 0
# dmax = 60

cmint = -3.
cmaxt = 3.
cmine = -0.2
cmaxe = 0.2
cminr = -0.5
cmaxr = 0.5
dmin = 0
dmax = 50

cmap = cm.balance

style.use('classic')
mpl.rcParams.update({'font.size': 8})
gs = gridspec.GridSpec(5, 10)

fig = plt.figure()
ax1 = plt.subplot(gs[0, 0:7])
ax2 = plt.subplot(gs[1:3, 0:7])
ax3 = plt.subplot(gs[3:5, 0:7])
ax4 = plt.subplot(gs[1:3, 8:11])
ax5 = plt.subplot(gs[3:5, 8:11])

ax1.grid('on')
ax4.grid('on')
ax5.grid('on')

ax1.set_xlim(yearday[0], yearday[-1])
ax2.set_xlim(yearday[0], yearday[-1])
ax3.set_xlim(yearday[0], yearday[-1])

ax1.set_ylim(cmint, cmaxt)
ax2.set_ylim(dmin, dmax)
ax3.set_ylim(dmin, dmax)
ax4.set_ylim(dmin, dmax)
ax5.set_ylim(dmin, dmax)

ax4.set_xlim(cmine, cmaxe)
ax5.set_xlim(cmine, cmaxe)

ax4.set_xticks([cmine, 0, cmaxe])
ax5.set_xticks([cmine, 0, cmaxe])

ax4.yaxis.tick_right()
ax5.yaxis.tick_right()
ax4.yaxis.set_label_position('right')
ax5.yaxis.set_label_position('right')

ax2.invert_yaxis()
ax3.invert_yaxis()
ax4.invert_yaxis()
ax5.invert_yaxis()

ax1.set_xticklabels([''])
ax2.set_xticklabels([''])
ax4.set_xticklabels([''])

ax1.set_ylabel(r'$U_{tide}$ [m$\cdot$s$^{-1}$]')
ax2.set_ylabel(r'Depth [m]')
ax3.set_xlabel(r'Yearday')
ax3.set_ylabel(r'Depth [m]')
ax4.set_ylabel(r'Depth [m]')
ax5.set_xlabel(r'$U_{exchange}$ [m$\cdot$s$^{-1}$]')
ax5.set_ylabel(r'Depth [m]')

ax1.plot(yearday, utide, 'b', label='U')
ax1.plot(yearday, vtide, 'r', label='V')
ax1.plot(yearday, np.zeros(yearday.shape), '--k')
ax1.plot(yearday2+0.5, vtidemax, '--k')
ax1.plot(yearday2+0.5, vtidemin, '--k')
ax1.legend(loc=(1.1, 0.2))

# pcm = ax2.contourf(yearday, z, u, np.linspace(cmin, cmax, 101),
#                    cmap=cmap, vmin=cmin, vmax=cmax, extend='both')
# pcm = ax3.contourf(yearday, z, v, np.linspace(cmin, cmax, 101),
#                    cmap=cmap, vmin=cmin, vmax=cmax, extend='both')
pcm = ax2.contourf(yearday, z, ures, np.linspace(cminr, cmaxr, 101),
                   cmap=cmap, vmin=cminr, vmax=cmaxr, extend='both')
pcm = ax3.contourf(yearday, z, vres, np.linspace(cminr, cmaxr, 101),
                   cmap=cmap, vmin=cminr, vmax=cmaxr, extend='both')
ax2.text(yearday[0]+1, 5, 'U')
ax3.text(yearday[0]+1, 5, 'V')

ax4.plot(utbar, z, 'b')
ax5.plot(vtbar, z, 'r')
ax4.plot([0, 0], [dmin, dmax], '--k')
ax5.plot([0, 0], [dmin, dmax], '--k')

# add colorbar axis handle
cbar_ax = fig.add_axes([0.67, 0.1, 0.01, 0.63])
cb = fig.colorbar(pcm, cax=cbar_ax, ticks=np.linspace(cminr, cmaxr, 11))
cbar_ax.set_ylabel(r'$U_{res}$ [m$\cdot$s$^{-1}$]')

plt.savefig('../figs/adcp_' + info['stn'] + '.png', dpi=600)
plt.close()
