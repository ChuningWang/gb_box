from datetime import datetime
import numpy as np
from scipy import interpolate
import netCDF4 as nc
from matplotlib import path
from ocean_toolbox import ctd
from gb_box import boxPrep

import matplotlib.pyplot as plt

box_info = {'case': 'Box',
            'box_method': 1,
            'clim': True,
            't0': nc.date2num(datetime(1990, 01, 01), 'days since 1900-01-01'),
            't1': nc.date2num(datetime(2010, 01, 01), 'days since 1900-01-01'),
            'grid_file_name': '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_lr_grd.nc',
            'river_file_name': '/Users/CnWang/git/gb_box/data/gb_box_rivers.nc',
            'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
            'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc'}

box1 = boxPrep.Box(box_info)
box1()

info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
        'filter': 'no',
       }

gb_ctd = ctd.ctd(info)
gb_ctd()

H_SI = 10
H_ID = 50

# -----------------------------------------------------------------------
# process river data
timer = box1.BoxRivers.climatology['time']
F0raw = box1.BoxRivers.climatology['river0']
F1raw = box1.BoxRivers.climatology['river1']
F2raw = box1.BoxRivers.climatology['river2']

# monthly average
F0 = np.zeros(12)
F1 = np.zeros(12)
F2 = np.zeros(12)
pytimer = nc.num2date(timer, 'days since 2000-01-01')
month = np.array([i.month for i in pytimer])
for i in range(1, 13):
    F0[i-1] = F0raw[month == i].mean()
    F1[i-1] = F1raw[month == i].mean()
    F2[i-1] = F2raw[month == i].mean()

# convert runoff unit from m3s-1 to m3day-1
F0 = F0*box1.consts['s2d']
F1 = F1*box1.consts['s2d']
F2 = F2*box1.consts['s2d']

# -----------------------------------------------------------------------
# process salinity data
time = gb_ctd.climatology['time']

lat_ctd = gb_ctd.data_info['lat'][gb_ctd.info['clim_station']]
lon_ctd = gb_ctd.data_info['lon'][gb_ctd.info['clim_station']]

coords = []
for i in range(len(lat_ctd)):
    coords.append((lon_ctd[i], lat_ctd[i]))

box_list = box1.boxes.keys()
path_points = {}
stn = {}
for box in box_list:
    path_points[box] = path.Path(box1.boxes[box])
    stn[box] = path_points[box].contains_points(coords)

salt0 = gb_ctd.climatology['salt'][:, :, stn['box0']]
salt0 = np.nanmean(salt0, axis=2)
salt1 = gb_ctd.climatology['salt'][:, :, stn['box1']]
salt1 = np.nanmean(salt1, axis=2)

Su0 = np.nanmean(salt0[:H_SI, :], axis=0)
Si0 = np.nanmean(salt0[H_SI:H_SI+H_ID, :], axis=0)
Sd0 = np.nanmean(salt0[H_SI+H_ID:, :], axis=0)

Su1 = np.nanmean(salt1[:H_SI, :], axis=0)
Si1 = np.nanmean(salt1[H_SI:H_SI+H_ID, :], axis=0)

# load box2 data
f = open('data/gb_icy_strait.txt', 'r')
text = f.readlines()
f.close()
text = [x.strip().split(',') for x in text]
time2 = np.array(text[0][1:]).astype(float)
Su2 = np.array(text[1][1:]).astype(float)
Si2 = np.array(text[2][1:]).astype(float)

# spline interpolate
spl = interpolate.splrep(time2, Su2)
Su2 = interpolate.splev(time, spl)
spl = interpolate.splrep(time2, Si2)
Si2 = interpolate.splev(time, spl)

# -----------------------------------------------------------------------
# calcualte d/dt
time = np.concatenate((time, [time[0] + 365]))
Su0 = np.concatenate((Su0, [Su0[0]]))
Si0 = np.concatenate((Si0, [Si0[0]]))
Sd0 = np.concatenate((Sd0, [Sd0[0]]))
Su1 = np.concatenate((Su1, [Su1[0]]))
Si1 = np.concatenate((Si1, [Si1[0]]))
Su2 = np.concatenate((Su2, [Su2[0]]))
Si2 = np.concatenate((Si2, [Si2[0]]))

F0 = np.concatenate((F0, [F0[0]]))
F1 = np.concatenate((F1, [F1[0]]))
F2 = np.concatenate((F2, [F2[0]]))

dt = np.diff(time)
dSu0dt = np.diff(Su0)/dt
dSi0dt = np.diff(Si0)/dt
dSd0dt = np.diff(Sd0)/dt
dSu1dt = np.diff(Su1)/dt
dSi1dt = np.diff(Si1)/dt
dSu2dt = np.diff(Su2)/dt
dSi2dt = np.diff(Si2)/dt

time = 0.5*(time[:-1] + time[1:])
Su0 = 0.5*(Su0[:-1] + Su0[1:])
Si0 = 0.5*(Si0[:-1] + Si0[1:])
Sd0 = 0.5*(Sd0[:-1] + Sd0[1:])
Su1 = 0.5*(Su1[:-1] + Su1[1:])
Si1 = 0.5*(Si1[:-1] + Si1[1:])
Su2 = 0.5*(Su2[:-1] + Su2[1:])
Si2 = 0.5*(Si2[:-1] + Si2[1:])

F0 = 0.5*(F0[:-1] + F0[1:])
F1 = 0.5*(F1[:-1] + F1[1:])
F2 = 0.5*(F2[:-1] + F2[1:])

# -----------------------------------------------------------------------
# construct the inverse matrix
month1 = 5
month2 = 10

Ainv = np.zeros(((month2-month1)*5, 5))
Binv = np.zeros((month2-month1)*5)

for i0, i in enumerate(range(5, 10)):
    Ainv[i0*5:(i0+1)*5, :] = np.array([
        [(Si0[i]-Su0[i])*(Su1[i]-Su0[i]),   0,                                  Si0[i]-Su0[i],      0,                  0],
        [(Si1[i]-Si0[i])*(Su1[i]-Su0[i]),   0,                                  -(Si0[i]-Su0[i]),   0,                  Sd0[i]-Si0[i]],
        [0,                                 0,                                  0,                  0,                  -(Sd0[i]-Si0[i])],
        [(Su0[i]-Su1[i])*(Su1[i]-Su0[i]),   (Si1[i]-Su1[i])*(Su2[i]-Su1[i]),    0,                  Si1[i]-Su1[i],      0],
        [(Su1[i]-Si1[i])*(Su1[i]-Su0[i]),   (Si2[i]-Si1[i])*(Su2[i]-Su1[i]),    0,                  -(Si1[i]-Su1[i]),   0]
        ])

    Binv[i0*5:(i0+1)*5] = np.array(
            [box1.volumes['Vs0']*dSu0dt[i] + F0[i]*Si0[i],
             box1.volumes['Vi0']*dSi0dt[i] - F0[i]*Si0[i] + F0[i]*Si1[i],
             box1.volumes['Vd0']*dSd0dt[i],
             box1.volumes['Vs1']*dSu1dt[i] + (F0[i] + F1[i])*Si1[i] - F0[i]*Su1[i],
             box1.volumes['Vi1']*dSi1dt[i] - (F0[i] + F1[i])*Si1[1] + (F0[i] + F1[i])*Si2[i]
                 - F0[i]*Si1[i] + F0[i]*Su1[i]])

[slv, _, _, _] = np.linalg.lstsq(Ainv, Binv)


# -----------------------------------------------------------------------
# save the data to txt file
fh = open('./data/gb_salt_clim' + str(box_info['box_method']) + '.txt', 'w')
fh.write('time, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(time))
fh.write('Su0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Su0))
fh.write('Si0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Si0))
fh.write('Sd0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Sd0))
fh.write('Su1, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Su1))
fh.write('Si1, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Si1))
fh.write('Su2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Su2))
fh.write('Si2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(Si2))
fh.write('dSu0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSu0dt))
fh.write('dSi0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSi0dt))
fh.write('dSd0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSd0dt))
fh.write('dSu1, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSu1dt))
fh.write('dSi1, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSi1dt))
fh.write('dSu2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSu2dt))
fh.write('dSi2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n' % tuple(dSi2dt))
fh.close()

# -----------------------------------------------------------------------
# construct the inverse matrix

# slv = []
# 
# for i in range(12):
#     Ainv = np.array([[Si0[i]-Su0[i],    0,              Si0[i]-Su0[i],      0,                  0],
#                      [Si1[i]-Si0[i],    0,              -(Si0[i]-Su0[i]),   0,                  Sd0[i]-Si0[i]],
#                      [0,                0,              0,                  0,                  -(Sd0[i]-Si0[i])],
#                      [Su0[i]-Su1[i],    Si1[i]-Su1[i],  0,                  Si1[i]-Su1[i],      0],
#                      [Su1[i]-Si1[i],    Si2[i]-Si1[i],  0,                  -(Si1[i]-Su1[i]),   0]
#                     ])
# 
#     A2inv = np.array([
#                      ])
# 
#     Binv = np.array([box1.volumes['Vs0']*dSu0dt[i] + F0[i]*Si0[i],
#                      box1.volumes['Vi0']*dSi0dt[i] - F0[i]*Si0[i] + F0[i]*Si1[i],
#                      box1.volumes['Vd0']*dSd0dt[i],
#                      box1.volumes['Vs1']*dSu1dt[i] + (F0[i] + F1[i])*Si1[i] - F0[i]*Su1[i],
#                      box1.volumes['Vi1']*dSi1dt[i] - (F0[i] + F1[i])*Si1[1] + (F0[i] + F1[i])*Si2[i]
#                          - F0[i]*Si1[i] + F0[i]*Su1[i],
#                     ])
# 
#     slv.append(np.linalg.solve(Ainv, Binv))
# 
# slv = np.array(slv)

# -----------------------------------------------------------------------
# # make a simple plot
# plt.plot(Su0)
# plt.plot(Si0)
# plt.plot(Sd0)
# plt.plot(Su1)
# plt.plot(Si1)
# plt.plot(Su2)
# plt.plot(Si2)
# plt.legend(['Su0', 'Si0', 'Sd0', 'Su1', 'Si1', 'Su2', 'Si2'])

# plt.plot(F0)
# plt.plot(F1)
# plt.plot(F2)
# plt.legend(['F0', 'F1', 'F2'])

# plt.plot(box1.BoxRivers.climatology['river0'])
# plt.plot(box1.BoxRivers.climatology['river1'])
# plt.plot(box1.BoxRivers.climatology['river2'])
# plt.legend(['river0', 'river1', 'river2'])
