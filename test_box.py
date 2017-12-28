""" test gb_box. """

from datetime import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from gb_box import boxPrep
from gb_box import boxODE

INFO = {'case': 'Box',
        'box_method': 1,
        'clim': True,
        'compare': True,
        't0': nc.date2num(datetime(2008, 01, 01), 'days since 1900-01-01'),
        't1': nc.date2num(datetime(2009, 01, 01), 'days since 1900-01-01'),
        'grid_file_name': '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_lr_grd.nc',
        'river_file_name': '/Users/CnWang/git/gb_box/data/gb_box_rivers.nc',
        'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
        'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc'}

box = boxPrep.Box(INFO)
box()

boxSolver = boxODE.BoxSolver()
boxSolver(box)

# make plots
fig, axarr = plt.subplots(4, sharex=True)
fig.subplots_adjust(hspace=0.05)
axarr[0].set_ylabel(r'Runoff [m$^3$s$^{-1}$]')
axarr[0].set_ylim(0, 1500)
axarr[0].set_yticks([0, 500, 1000, 1500])
axarr[2].set_ylabel(r'Salinity [PSU]')
axarr[3].set_xlim(0, 366)
axarr[3].set_xlabel('Yearday')
axarr[1].set_ylim(23, 33)
axarr[1].set_yticks([24, 26, 28, 30, 32])
axarr[2].set_ylim(23, 33)
axarr[2].set_yticks([24, 26, 28, 30, 32])
axarr[3].set_ylim(23, 33)
axarr[3].set_yticks([24, 26, 28, 30, 32])

axarr[0].plot(box.time, box.BoxRivers.data['river0'], 'b', label='River 0')
axarr[0].plot(box.time, box.BoxRivers.data['river1'], 'r', label='River 1')
axarr[0].plot(box.time, box.BoxRivers.data['river2'], 'k', label='River 2')
axarr[0].legend()

axarr[1].plot(box.measurements['time'], box.measurements['Su0'],
              '--.b', linewidth=.3)
axarr[1].plot(box.measurements['time'], box.measurements['Si0'],
              '--.r', linewidth=.3)
axarr[1].plot(box.measurements['time'], box.measurements['Sd0'],
              '--.k', linewidth=.3)
axarr[1].legend((r'$S_{u0}$', r'$S_{i0}$', r'$S_{d0}$'), loc=3)

axarr[2].plot(box.measurements['time'], box.measurements['Su1'],
              '--.b', linewidth=.3)
axarr[2].plot(box.measurements['time'], box.measurements['Si1'],
              '--.r', linewidth=.3)
axarr[2].legend((r'$S_{u1}$', r'$S_{i1}$'), loc=3)

axarr[3].plot(box.measurements['time'], box.measurements['Su2'],
              '--.b', linewidth=.3)
axarr[3].plot(box.measurements['time'], box.measurements['Si2'],
              '--.r', linewidth=.3)
axarr[3].legend((r'$S_{u2}$', r'$S_{i2}$'), loc=3)

axarr[1].plot(box.time, boxSolver.solu[0, :]*box.consts['s0'], '--b')
axarr[1].plot(box.time, boxSolver.solu[1, :]*box.consts['s0'], '--r')
axarr[1].plot(box.time, boxSolver.solu[2, :]*box.consts['s0'], '--k')

axarr[2].plot(box.time, boxSolver.solu[3, :]*box.consts['s0'], '--b')
axarr[2].plot(box.time, boxSolver.solu[4, :]*box.consts['s0'], '--r')

axarr[3].plot(box.time, boxSolver.solu[6, :]*box.consts['s0'], '--b')
axarr[3].plot(box.time, boxSolver.solu[7, :]*box.consts['s0'], '--r')

plt.savefig('salt.png', dpi=300)
plt.close()
# plt.show()
