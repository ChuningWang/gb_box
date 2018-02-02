import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import style

from cmocean import cm
from ocean_toolbox import ctd

info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd_clim_all.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': range(25),
        'clim_deep_interp': 'no',
        'filter': 'no'}
gb_ctd = ctd.ctd(info)
gb_ctd()
stn_list = [21, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 24]
gb_ctd.get_trans_clim(['salt', 'temp'], stn_list)

# ------------ make plots ----------------------------------------------------
depth0 = 50
depth1 = 450
dtrans = 150

salt1 = gb_ctd.trans['salt'][:depth1, 0, :]
salt2 = gb_ctd.trans['salt'][:depth1, 6, :]
temp1 = gb_ctd.trans['temp'][:depth1, 0, :]
temp2 = gb_ctd.trans['temp'][:depth1, 6, :]

salt1 = np.vstack((salt1, np.nan*np.zeros((1, len(stn_list)))))
salt2 = np.vstack((salt2, np.nan*np.zeros((1, len(stn_list)))))
temp1 = np.vstack((temp1, np.nan*np.zeros((1, len(stn_list)))))
temp2 = np.vstack((temp2, np.nan*np.zeros((1, len(stn_list)))))

depth = gb_ctd.trans['fathometer_depth']
depth = (depth-depth0)*(depth1-dtrans)/(depth1-depth0) + dtrans

z1 = np.linspace(0, dtrans, depth0+1)
z2 = np.linspace(dtrans, depth1, depth1-depth0+1)
znew = np.concatenate((z1[:-1], z2))
yticknew = np.concatenate((np.linspace(0, dtrans, 6)[:-1], np.linspace(dtrans, depth1, 9)))
yticktrue = np.concatenate((np.linspace(0, depth0, 6)[:-1], np.linspace(depth0, depth1, 9)))
yticktrue = yticktrue.astype(int)
yticktrue = map(str, yticktrue)

style.use('classic')
mpl.rcParams.update({'font.size': 8})

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.1, wspace=0.25)

ax[0, 0].set_xlim(gb_ctd.trans['dis'][0], gb_ctd.trans['dis'][-1])
ax[0, 0].set_ylim(0, 450)
ax[0, 0].set_yticks(yticknew)
ax[0, 0].set_yticklabels(yticktrue)
ax[0, 0].invert_yaxis()

ax[1, 0].set_xlabel('Distance [km]')
ax[1, 1].set_xlabel('Distance [km]')
ax[0, 0].set_ylabel('Depth [m]')
ax[1, 0].set_ylabel('Depth [m]')

# contourf variable
ctf0 = ax[0, 0].contourf(gb_ctd.trans['dis'], znew,
                         salt1, np.linspace(12, 32, 101),
                         vmin=12, vmax=32, cmap=cm.haline, extend='both')

# plot bathymetry
ax[0, 0].fill_between(gb_ctd.trans['dis'], depth,
                      450, facecolor='lightgrey')
ax[0, 0].plot(gb_ctd.trans['dis'], 150*np.ones(len(stn_list)), '--k', lw=1.5)

# contourf variable
ctf0 = ax[1, 0].contourf(gb_ctd.trans['dis'], znew,
                         salt2, np.linspace(12, 32, 101),
                         vmin=12, vmax=32, cmap=cm.haline, extend='both')

# plot bathymetry
ax[1, 0].fill_between(gb_ctd.trans['dis'], depth,
                      450, facecolor='lightgrey')
ax[1, 0].plot(gb_ctd.trans['dis'], 150*np.ones(len(stn_list)), '--k', lw=1.5)

# add colorbar axis handle
cbar_ax = fig.add_axes([0.48, 0.1, 0.01, 0.8])
cb = fig.colorbar(ctf0, cax=cbar_ax, ticks=np.linspace(12, 32, 11))
cbar_ax.set_ylabel('Salinity [PSU]')

ax[0, 0].text(5, 420, 'a) Jan', fontsize=10)
ax[1, 0].text(5, 420, 'b) Jul', fontsize=10)

# contourf variable
ctf0 = ax[0, 1].contourf(gb_ctd.trans['dis'], znew,
                         temp1, np.linspace(2, 12, 101),
                         vmin=2, vmax=12, cmap=cm.thermal, extend='both')

# plot bathymetry
ax[0, 1].fill_between(gb_ctd.trans['dis'], depth,
                      450, facecolor='lightgrey')
ax[0, 1].plot(gb_ctd.trans['dis'], 150*np.ones(len(stn_list)), '--k', lw=1.5)

# contourf variable
ctf0 = ax[1, 1].contourf(gb_ctd.trans['dis'], znew,
                         temp2, np.linspace(2, 12, 101),
                         vmin=2, vmax=12, cmap=cm.thermal, extend='both')

# plot bathymetry
ax[1, 1].fill_between(gb_ctd.trans['dis'], depth,
                      450, facecolor='lightgrey')
ax[1, 1].plot(gb_ctd.trans['dis'], 150*np.ones(len(stn_list)), '--k', lw=1.5)

# add colorbar axis handle
cbar_ax = fig.add_axes([0.91, 0.1, 0.01, 0.8])
cb = fig.colorbar(ctf0, cax=cbar_ax, ticks=np.linspace(2, 12, 11))
cbar_ax.set_ylabel('Temperature [$^{\circ}$C]')

ax[0, 1].text(5, 420, 'c) Jan', fontsize=10)
ax[1, 1].text(5, 420, 'd) Jul', fontsize=10)

# plot station location
for i in [0, 1]:
    for j in [0, 1]:
        ax[i, j].plot(gb_ctd.trans['dis'], np.zeros(gb_ctd.trans['dis'].shape),
                      'vk', markersize=5)
        for k, stni in enumerate(stn_list):
            ax[i, j].text(gb_ctd.trans['dis'][k], -20, '%02d' % stni,
                          horizontalalignment='center', fontsize=7, rotation=60)
plt.grid('on', linewidth=0.1)
ax[0, 0].text(40, -40, 'Station Number')
ax[0, 1].text(40, -40, 'Station Number')

plt.savefig('figs/clim_all2.png', dpi=600)
plt.close()


# var = 'temp'
# cmin = 2
# cmid = 5
# cmax = 12
# cmap = cm.thermal
# # longname = r'Salinity [PSU]'
# longname = r'Temperature [$^{\circ}$C]'
# data = gb_ctd.trans[var]
# data = np.ma.masked_invalid(data)

# clevs = np.concatenate((np.arange(cmin, cmid, 5), np.linspace(cmid, cmax, 11)))
# cb_ticks = [5, 10, 15, 20, 25, 30, 30.5, 31, 31.5, 32]
# clevs = np.linspace(cmin, cmax, 101)
# cb_ticks = clevs[::10]

# fig, ax = plt.subplots(4, 2, sharex=True)
# fig.subplots_adjust(hspace=0.05, wspace=0.2)
# ax[0, 0].set_xlim(gb_ctd.trans['dis'][0], gb_ctd.trans['dis'][-1])
# ax[0, 0].set_ylim(0, depth0)
# ax[0, 0].set_yticks(range(0, depth0, 10))
# ax[0, 0].invert_yaxis()
# ax[1, 0].set_ylim(depth0, depth1)
# ax[1, 0].set_yticks(range(depth0, depth1, 100))
# ax[1, 0].invert_yaxis()
# 
# ax[2, 0].set_xlim(gb_ctd.trans['dis'][0], gb_ctd.trans['dis'][-1])
# ax[2, 0].set_ylim(0, depth0)
# ax[2, 0].set_yticks(range(0, depth0, 10))
# ax[2, 0].invert_yaxis()
# ax[3, 0].set_ylim(depth0, depth1)
# ax[3, 0].set_yticks(range(depth0, depth1, 100))
# ax[3, 0].invert_yaxis()
# 
# ax[0, 1].set_xlim(gb_ctd.trans['dis'][0], gb_ctd.trans['dis'][-1])
# ax[0, 1].set_ylim(0, depth0)
# ax[0, 1].set_yticks(range(0, depth0, 10))
# ax[0, 1].invert_yaxis()
# ax[1, 1].set_ylim(depth0, depth1)
# ax[1, 1].set_yticks(range(depth0, depth1, 100))
# ax[1, 1].invert_yaxis()
# 
# ax[2, 1].set_xlim(gb_ctd.trans['dis'][0], gb_ctd.trans['dis'][-1])
# ax[2, 1].set_ylim(0, depth0)
# ax[2, 1].set_yticks(range(0, depth0, 10))
# ax[2, 1].invert_yaxis()
# ax[3, 1].set_ylim(depth0, depth1)
# ax[3, 1].set_yticks(range(depth0, depth1, 100))
# ax[3, 1].invert_yaxis()
# 
# ax[3, 0].set_xlabel('Distance [km]')
# ax[3, 0].set_ylabel('Depth [m]')
# 
# ax[0, 1].set_yticklabels({''})
# ax[1, 1].set_yticklabels({''})
# ax[2, 1].set_yticklabels({''})
# ax[3, 1].set_yticklabels({''})
# 
# # contourf variable
# ctf0 = ax[0, 0].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['salt'][:, 0, :], np.linspace(15, 30, 101),
#                          vmin=15, vmax=30, cmap=cm.haline, extend='both')
# ctf1 = ax[1, 0].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['salt'][:, 0, :], np.linspace(15, 30, 101),
#                          vmin=15, vmax=30, cmap=cm.haline, extend='both')
# 
# # plot bathymetry
# ax[0, 0].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# ax[1, 0].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# 
# # contourf variable
# ctf0 = ax[2, 0].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['salt'][:, 6, :], np.linspace(15, 30, 101),
#                          vmin=15, vmax=30, cmap=cm.haline, extend='both')
# ctf1 = ax[3, 0].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['salt'][:, 6, :], np.linspace(15, 30, 101),
#                          vmin=15, vmax=30, cmap=cm.haline, extend='both')
# 
# # plot bathymetry
# ax[2, 0].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# ax[3, 0].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# 
# # add colorbar axis handle
# cbar_ax = fig.add_axes([0.48, 0.1, 0.01, 0.8])
# cb = fig.colorbar(ctf1, cax=cbar_ax, ticks=np.linspace(15, 30, 16))
# cbar_ax.set_ylabel('Salinity [PSU]')
# 
# # contourf variable
# ctf0 = ax[0, 1].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['temp'][:, 0, :], np.linspace(2, 12, 101),
#                          vmin=2, vmax=12, cmap=cm.thermal, extend='both')
# ctf1 = ax[1, 1].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['temp'][:, 0, :], np.linspace(2, 12, 101),
#                          vmin=2, vmax=12, cmap=cm.thermal, extend='both')
# 
# # plot bathymetry
# ax[0, 1].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# ax[1, 1].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# 
# # contourf variable
# ctf0 = ax[2, 1].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['temp'][:, 6, :], np.linspace(2, 12, 101),
#                          vmin=2, vmax=12, cmap=cm.thermal, extend='both')
# ctf1 = ax[3, 1].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                          gb_ctd.trans['temp'][:, 6, :], np.linspace(2, 12, 101),
#                          vmin=2, vmax=12, cmap=cm.thermal, extend='both')
# 
# # plot bathymetry
# ax[2, 1].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# ax[3, 1].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                       depth1, facecolor='lightgrey')
# 
# # add colorbar axis handle
# cbar_ax = fig.add_axes([0.9025, 0.1, 0.01, 0.8])
# cb = fig.colorbar(ctf1, cax=cbar_ax, ticks=np.linspace(2, 12, 11))
# cbar_ax.set_ylabel(r'Temperature [$^{\circ}$C]')
# 
# ax[1, 0].text(5, 380, 'Jan')
# ax[1, 1].text(5, 380, 'Jan')
# ax[3, 0].text(5, 380, 'Jul')
# ax[3, 1].text(5, 380, 'Jul')
# 
# plt.savefig('figs/clim_all.png', dpi=600)
# plt.close()

# for mm in range(12):
# 
#     # set the axises
#     fig, ax = plt.subplots(2, sharex=True)
#     fig.subplots_adjust(hspace=0.05)
#     ax[0].set_xlim(gb_ctd.trans['dis'][0], gb_ctd.trans['dis'][-1])
#     ax[0].set_ylim(0, depth0)
#     ax[0].set_yticks(range(0, depth0, 10))
#     ax[0].invert_yaxis()
#     ax[1].set_ylim(depth0, depth1)
#     ax[1].set_yticks(range(depth0, depth1, 100))
#     ax[1].invert_yaxis()
#     # labels
#     ax[1].set_xlabel('Distance [km]')
#     ax[0].set_ylabel('Depth [m]')
# 
#     # set axis position
#     pos = ax[0].get_position()
#     pos2 = [pos.x0-0.04, pos.y0,  pos.width, pos.height]
#     ax[0].set_position(pos2)
#     pos = ax[1].get_position()
#     pos2 = [pos.x0-0.04, pos.y0,  pos.width, pos.height]
#     ax[1].set_position(pos2)
# 
#     # contourf variable
#     ctf0 = ax[0].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                           data[:, mm, :], clevs,
#                           vmin=cmin, vmax=cmax, cmap=cmap, extend='both')
#     ctf1 = ax[1].contourf(gb_ctd.trans['dis'], gb_ctd.trans['z'],
#                           data[:, mm, :], clevs,
#                           vmin=cmin, vmax=cmax, cmap=cmap, extend='both')
# 
#     # plot bathymetry
#     ax[0].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                        depth1, facecolor='lightgrey')
#     ax[1].fill_between(gb_ctd.trans['dis'], gb_ctd.trans['fathometer_depth'],
#                        depth1, facecolor='lightgrey')
# 
#     # add colorbar axis handle
#     cbar_ax = fig.add_axes([0.88, 0.1, 0.02, 0.8])
#     cb = fig.colorbar(ctf1, cax=cbar_ax, ticks=cb_ticks)
#     cbar_ax.set_ylabel(longname)
# 
#     # plot station location
#     ax[0].text(45, -8, 'Station Number')
#     ax[0].plot(gb_ctd.trans['dis'], np.zeros(gb_ctd.trans['dis'].shape),
#                'vk', markersize=5)
#     plt.grid('on', linewidth=0.1)
#     for i, stni in enumerate(stn_list):
#         ax[0].text(gb_ctd.trans['dis'][i], -2, '%02d' % stni,
#                    horizontalalignment='center')
#     # ttl = var + '_' + nc.num2date(self.trans['time'], 'days since 1900-01-01').strftime('%Y-%m-%d')
#     # ax[0].set_title(ttl)
# 
#     plt.savefig('figs/' + var + '_clim_%02d.png' % (mm+1), dpi=600)
#     plt.close()
