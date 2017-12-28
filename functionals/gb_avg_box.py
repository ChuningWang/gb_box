""" calculate box averaged salinity climatology (monthly). """

from datetime import datetime
import os
import sys

import numpy as np
from scipy import interpolate
from matplotlib import path
import netCDF4 as nc

sys.path.insert(0, os.path.abspath(".."))

from ocean_toolbox import ctd
from gb_box import boxPrep

box_info = {'case': 'Box2',
            'box_method': 2,
            'clim': True,
            't0': nc.date2num(datetime(1990, 01, 01), 'days since 1900-01-01'),
            't1': nc.date2num(datetime(2010, 01, 01), 'days since 1900-01-01'),
            'grid_file_name': '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_lr_grd.nc',
            'river_file_name': '/Users/CnWang/git/gb_box/data/gb_box_rivers2.nc',
            'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
            'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc'}

box = boxPrep.Box(box_info)
box()

info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
        'filter': 'no'}

gb_ctd = ctd.ctd(info)
gb_ctd()

H_SI = 10
H_ID = 50

# -----------------------------------------------------------------------
# process salinity data
time = gb_ctd.climatology['time']

lat_ctd = gb_ctd.data_info['lat'][gb_ctd.info['clim_station']]
lon_ctd = gb_ctd.data_info['lon'][gb_ctd.info['clim_station']]

coords = []
for i in range(len(lat_ctd)):
    coords.append((lon_ctd[i], lat_ctd[i]))

path_points = {}
stn = {}
for boxi in box.boxes_list:
    path_points[boxi] = path.Path(box.boxes[boxi])
    stn[boxi] = path_points[boxi].contains_points(coords)

Su = {}
Si = {}
Sd = {}
for boxi in box.boxes_list:
    salt = gb_ctd.climatology['salt'][:, :, stn[boxi]]
    salt = np.nanmean(salt, axis=2)
    Su[boxi] = np.nanmean(salt[:H_SI, :], axis=0)
    Si[boxi] = np.nanmean(salt[H_SI:H_SI+H_ID, :], axis=0)
    Sd[boxi] = np.nanmean(salt[H_SI+H_ID:, :], axis=0)

# -----------------------------------------------------------------------
# load box4 data
fin = open('../data/gb_icy_strait.txt', 'r')
text = fin.readlines()
fin.close()
text = [x.strip().split(',') for x in text]
time2 = np.array(text[0][1:]).astype(float)
Su4 = np.array(text[1][1:]).astype(float)
Si4 = np.array(text[2][1:]).astype(float)

# spline interpolate
spl = interpolate.splrep(time2, Su4)
Su['box4'] = interpolate.splev(time, spl)
spl = interpolate.splrep(time2, Si4)
Si['box4'] = interpolate.splev(time, spl)

Sd['box4'][:] = np.NaN

# -----------------------------------------------------------------------
# combine some boxes, mask some invalid boxes
Sdc = (Sd['box0']*box.volumes['Vd0'] + Sd['box1']*box.volumes['Vd1'] + Sd['box2']*box.volumes['Vd2']) / \
    (box.volumes['Vd0'] + box.volumes['Vd1'] + box.volumes['Vd2'])

Sd['box2'] = Sdc

Sd['box0'][:] = np.NaN
Sd['box1'][:] = np.NaN
Sd['box3'][:] = np.NaN

# -----------------------------------------------------------------------
# calcualte d/dt
# time = np.concatenate(([time[-1] - 365], time, [time[0] + 365]))
# time = 0.5*(time[:-1] + time[1:])

dSudt = {}
dSidt = {}
dSddt = {}
for boxi in box.boxes_list:

    dSudt[boxi] = np.zeros(Su[boxi].shape)
    dSidt[boxi] = np.zeros(Si[boxi].shape)
    dSddt[boxi] = np.zeros(Sd[boxi].shape)

    dSudt[boxi][1:-1] = (Su[boxi][2:]-Su[boxi][:-2])/(time[2:]-time[:-2])
    dSidt[boxi][1:-1] = (Si[boxi][2:]-Si[boxi][:-2])/(time[2:]-time[:-2])
    dSddt[boxi][1:-1] = (Sd[boxi][2:]-Sd[boxi][:-2])/(time[2:]-time[:-2])

    dSudt[boxi][0] = (Su[boxi][1]-Su[boxi][-1])/(time[1]-(time[-1]-365))
    dSidt[boxi][0] = (Si[boxi][1]-Si[boxi][-1])/(time[1]-(time[-1]-365))
    dSddt[boxi][0] = (Sd[boxi][1]-Sd[boxi][-1])/(time[1]-(time[-1]-365))

    dSudt[boxi][-1] = (Su[boxi][0]-Su[boxi][-2])/(time[0]-(time[-2]-365))
    dSidt[boxi][-1] = (Si[boxi][0]-Si[boxi][-2])/(time[0]-(time[-2]-365))
    dSddt[boxi][-1] = (Sd[boxi][0]-Sd[boxi][-2])/(time[0]-(time[-2]-365))

# # -----------------------------------------------------------------------
# # # save the data to txt file
fmt = ', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n'
fh = open('../data/gb_salt_clim_' + str(box_info['box_method']) + '.txt', 'w')
fh.write('time' + fmt % tuple(time))
for boxi in box.boxes_list:
    fh.write('Su' + boxi[-1] + fmt % tuple(Su[boxi]))
    fh.write('Si' + boxi[-1] + fmt % tuple(Si[boxi]))
    fh.write('Sd' + boxi[-1] + fmt % tuple(Sd[boxi]))
for boxi in box.boxes_list:
    fh.write('dSu' + boxi[-1] + fmt % tuple(dSudt[boxi]))
    fh.write('dSi' + boxi[-1] + fmt % tuple(dSidt[boxi]))
    fh.write('dSd' + boxi[-1] + fmt % tuple(dSddt[boxi]))
fh.close()
