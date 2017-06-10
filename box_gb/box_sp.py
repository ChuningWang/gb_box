# Read off-coast salinity climatology and calculate Sp

import numpy as np
from scipy.interpolate import interp1d

def box_sp(mt, lat_sp=57.875, lon_sp=-136.625, d_lim = [50, 100]):

    import netCDF4 as nc

    # lat_sp = 57.875
    # lon_sp = -136.625
    # d_lim = [20, 100]

    import matplotlib.dates as mdates
    from datetime import datetime
    import statsmodels.api as sm
    # import pdb

    # Load WOD salinity data
    fh = nc.Dataset('./data/s_wod.nc', 'r')
    mm = fh.variables['mm'][:]
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    d = fh.variables['d'][:]
    s = fh.variables['s'][:]

    # Select data at coordinate (lon_sp, lat_sp)
    msk1 = lat==lat_sp
    msk2 = lon==lon_sp
    msk3 = (d>=d_lim[0]) & (d<=d_lim[1])
    sp = s[:, msk3, msk1, msk2]
    sp = np.mean(sp, axis=1)

    yd = np.zeros(12)
    for i in range(yd.size):
        yd[i] = datetime(1, i+1, 15).timetuple().tm_yday

    yd = np.concatenate((yd-366, yd, yd+366))
    sp = np.concatenate((sp, sp, sp))

    yd_i = np.arange(-300, 300*2)
    sp_i = interp1d(yd,sp)(yd_i)

    # LOWESS filter
    lowess = sm.nonparametric.lowess(sp_i, yd_i, frac=0.2)
    # Interpolate onto mt
    # Since num2date cannot deal with mt<1, add 366 days (1 year) if mt<1
    mt2 = np.zeros(mt.size)
    mt2[:] = mt
    mt2[mt2<1] = mt2[mt2<1]+366
    pyt = np.array(mdates.num2date(mt2))
    yd_mt = np.array([pyt[i].timetuple().tm_yday for i in range(pyt.size)])
    sp = interp1d(lowess[:, 0], lowess[:, 1])(yd_mt)

    # pdb.set_trace()

    return sp

