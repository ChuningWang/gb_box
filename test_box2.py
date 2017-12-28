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
        # 'sl_rivers': 's',
        # 'river_raw_file_name': '/Users/CnWang/Documents/gb_roms/gb_discharge.nc',
        'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
        'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc'}

box = boxPrep.Box(INFO)
box()

boxSolver = boxODE.BoxSolver()
boxSolver(box)

# ---------- make plots ---------------------------------------------
fig, axarr = plt.subplots(5, sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.05, right = 0.87)
axarr2 = np.array([i.twinx() for i in axarr])
# axarr[0].set_ylim(0, 1500)
# axarr[0].set_yticks([0, 500, 1000, 1500])
axarr[2].set_ylabel(r'Salinity [PSU]')
axarr[-1].set_xlim(0, 366)
axarr[-1].set_xlabel('Yearday')
# axarr[-1].set_ylim(23, 34)
# axarr[-1].set_yticks([24, 27, 30, 33])

for i in axarr2:
    i.set_ylim(0, 1500)
    i.set_yticks([500, 1000, 1500])
    i.tick_params('y', colors='r')

axarr2[2].set_ylabel(r'Runoff [m$^3$s$^{-1}$]', color='r')

for i in range(5):
    axarr[i].plot(box.measurements['time'], box.measurements['Su'+str(i)],
                  '--.b', linewidth=.3, label='Su')
    axarr[i].plot(box.measurements['time'], box.measurements['Si'+str(i)],
                  '--.k', linewidth=.3, label='Si')
    axarr[i].plot(box.time, boxSolver.solution[3*i, :]*box.consts['s0'], '--b')
    axarr[i].plot(box.time, boxSolver.solution[3*i+1, :]*box.consts['s0'], '--k')
    axarr2[i].plot(box.time, box.BoxRivers.data['river' + str(i)],
                   'r', label='River 0')

axarr[0].plot(box.measurements['time'], box.measurements['Sd2'],
              '--.g', linewidth=.3, label='Sd')
axarr[1].plot(box.measurements['time'], box.measurements['Sd2'],
              '--.g', linewidth=.3, label='Sd')
axarr[2].plot(box.measurements['time'], box.measurements['Sd2'],
              '--.g', linewidth=.3, label='Sd')
axarr[0].plot(box.time, boxSolver.solution[8, :]*box.consts['s0'], '--g')
axarr[1].plot(box.time, boxSolver.solution[8, :]*box.consts['s0'], '--g')
axarr[2].plot(box.time, boxSolver.solution[8, :]*box.consts['s0'], '--g')

plt.savefig('salt2.png', dpi=300)
plt.close()

# # make plots
# fig, axarr = plt.subplots(4, sharex=True)
# fig.subplots_adjust(hspace=0.05)
# axarr[0].set_ylabel(r'Runoff [m$^3$s$^{-1}$]')
# axarr[0].set_ylim(0, 1500)
# axarr[0].set_yticks([0, 500, 1000, 1500])
# axarr[2].set_ylabel(r'Salinity [PSU]')
# axarr[3].set_xlim(0, 366)
# axarr[3].set_xlabel('Yearday')
# axarr[1].set_ylim(23, 33)
# axarr[1].set_yticks([24, 26, 28, 30, 32])
# axarr[2].set_ylim(23, 33)
# axarr[2].set_yticks([24, 26, 28, 30, 32])
# axarr[3].set_ylim(23, 33)
# axarr[3].set_yticks([24, 26, 28, 30, 32])
# 
# axarr[0].plot(box.time, box.BoxRivers.data['river0'], 'b', label='River 0')
# axarr[0].plot(box.time, box.BoxRivers.data['river1'], 'r', label='River 1')
# axarr[0].plot(box.time, box.BoxRivers.data['river2'], 'k', label='River 2')
# axarr[0].legend()
# 
# axarr[1].plot(box.measurements['time'], box.measurements['Su0'],
#               '--.b', linewidth=.3)
# axarr[1].plot(box.measurements['time'], box.measurements['Si0'],
#               '--.r', linewidth=.3)
# axarr[1].plot(box.measurements['time'], box.measurements['Sd0'],
#               '--.k', linewidth=.3)
# axarr[1].legend((r'$S_{u0}$', r'$S_{i0}$', r'$S_{d0}$'), loc=3)
# 
# axarr[2].plot(box.measurements['time'], box.measurements['Su1'],
#               '--.b', linewidth=.3)
# axarr[2].plot(box.measurements['time'], box.measurements['Si1'],
#               '--.r', linewidth=.3)
# axarr[2].legend((r'$S_{u1}$', r'$S_{i1}$'), loc=3)
# 
# axarr[3].plot(box.measurements['time'], box.measurements['Su2'],
#               '--.b', linewidth=.3)
# axarr[3].plot(box.measurements['time'], box.measurements['Si2'],
#               '--.r', linewidth=.3)
# axarr[3].legend((r'$S_{u2}$', r'$S_{i2}$'), loc=3)
# 
# axarr[1].plot(box.time, boxSolver.solu[0, :]*box.consts['s0'], '--b')
# axarr[1].plot(box.time, boxSolver.solu[1, :]*box.consts['s0'], '--r')
# axarr[1].plot(box.time, boxSolver.solu[2, :]*box.consts['s0'], '--k')
# 
# axarr[2].plot(box.time, boxSolver.solu[3, :]*box.consts['s0'], '--b')
# axarr[2].plot(box.time, boxSolver.solu[4, :]*box.consts['s0'], '--r')
# 
# axarr[3].plot(box.time, boxSolver.solu[6, :]*box.consts['s0'], '--b')
# axarr[3].plot(box.time, boxSolver.solu[7, :]*box.consts['s0'], '--r')
# 
# plt.savefig('salt.png', dpi=300)
# plt.close()
# # plt.show()
