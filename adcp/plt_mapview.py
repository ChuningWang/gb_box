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
# stn = 'SEA0850'
# stn = 'SEA1003'
# stn = 'SEA1004'

stn_list = ['SEA0830', 'SEA0831', 'SEA0832', 'SEA0839', 'SEA0840',
            'SEA0841', 'SEA0842', 'SEA0843', 'SEA0844', 'SEA0845',
            'SEA0847', 'SEA0850', 'SEA1001', 'SEA1002', 'SEA1003',
            'SEA1004', 'SEA1005', 'SEA1006', 'SEA1007', 'SEA1008',
            'SEA1009', 'SEA1010']

lon = []
lat = []
Ubar0 = []
Ut0 = []

# -------------- extract data -------------------------------
info = {'stn' : stn,
        'file_dir': out_dir + 'tef/NOAA_ADCP/',
        'sl': 'l',
        'Wp_hrs': 24}

crt = noaa_adcp.get_noaa_current(info)
crt()
time = crt.ctime
U = crt.uraw + 1j*crt.vraw
Ubar = U.mean(axis=0)
# time convertion
pytime = nc.num2date(time, 'days since 1900-01-01')
yearday = time - \
    nc.date2num(datetime(pytime[0].year, 1, 1), 'days since 1900-01-01') + \
    1
time2 = np.arange(np.floor(time[0]), np.floor(time[0]))
yearday2 = np.arange(np.floor(yearday[0]), np.floor(yearday[-1]))

lon.apend(crt.info['lon'])
lat.apend(crt.info['lat'])

# -------------- harmonic analysis --------------------------
tfit = ttide.t_tide(Ubar, dt=0.1,
                    stime=pytime[0], out_style=None, errcalc='wboot')
m2 = tfit['nameu'] == 'M2  '
ang = tfit['tidecon'][m2, 4][0]-90
Ut = tfit(pytime)

