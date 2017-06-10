# Box model for Glacer Bay
# The box model is developed based on the SoG box model. 
# Use to initiate the Box Model (generate IC and BC, etc.)
#
# Chuning Wang

# -----------------------------------------------------------------------------------------
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from scipy.interpolate import interp1d

import pdb

# -----------------------------------------------------------------------------------------
def get_area_volume(hu, hi, hd, pth='data/ARDEMv2.0.nc', boxMethod=1):
    # Calculate volume and area for each box
    from gb_toolbox.gb_discharge import get_avgbox
    from matplotlib import path
    from geopy.distance import vincenty
    import netCDF4 as nc

    print 'Calculating area and volume for each box...'
    depth = hu+hi+hd
    fh = nc.Dataset(pth, mode='r')
    lon = fh.variables['lon'][:]-360
    lat = fh.variables['lat'][:]
    z = fh.variables['z'][:]
    msk1 = (lat > 57.5) & (lat < 59.5)
    msk2 = (lon > -138) & (lon < -135)
    lat = lat[msk1]
    lon = lon[msk2]
    z = z[msk1, :][:, msk2]

    dis1 = np.zeros((lat.size-1, lon.size-1))
    dis2 = np.zeros((lat.size-1, lon.size-1))
    for i in range(lat.size-1):
        for j in range(lon.size-1):
            dis1[i, j] = vincenty((lat[i], lon[j]), (lat[i+1], lon[j])).meters
            dis2[i, j] = vincenty((lat[i], lon[j]), (lat[i], lon[j+1])).meters

    Area = dis1*dis2

    # Get points in boxes
    boxes = get_avgbox(boxMethod=boxMethod)
    topo_box = np.ones((lat.size, lon.size))*(-1)
    p0 = path.Path(boxes['box0'])
    p1 = path.Path(boxes['box1'])
    p2 = path.Path(boxes['box2'])
    p3 = path.Path(boxes['box3'])

    for i in range(lat.size):
        for j in range(lon.size):
            if p0.contains_points([(lon[j], lat[i])]):
                topo_box[i, j] = 0
            elif p1.contains_points([(lon[j], lat[i])]):
                topo_box[i, j] = 1
            elif p2.contains_points([(lon[j], lat[i])]):
                topo_box[i, j] = 2
            elif p3.contains_points([(lon[j], lat[i])]):
                topo_box[i, j] = 3

    # Remove the last column and row to match Area
    topo_box = topo_box[:-1, :-1]
    z = z[:-1, :-1]

    # Calculate basin areas for each layer (dz = 1m)
    A_layer = np.zeros((depth, 4))
    for i in range(depth):
        topo_box[z > -i] = -1

        # Calculate Layer Area
        A_layer[i, 0] = np.sum(Area[topo_box == 0])
        A_layer[i, 1] = np.sum(Area[topo_box == 1])
        A_layer[i, 2] = np.sum(Area[topo_box == 2])
        A_layer[i, 3] = np.sum(Area[topo_box == 3])

    # Get interface of Area between two layers
    As, Aui, Aid = np.zeros(4), np.zeros(4), np.zeros(4)
    As[:] = A_layer[0, :]
    Aui[:] = A_layer[hu, :]
    Aid[:] = A_layer[hu+hi, :]

    # Get volume of each box
    Vu, Vi, Vd = np.zeros(4), np.zeros(4), np.zeros(4)
    Vu[:] = np.sum(A_layer[:hu, :], axis=0)
    Vi[:] = np.sum(A_layer[hu:hu+hi, :], axis=0)
    Vd[:] = np.sum(A_layer[hu+hi:, :], axis=0)

    return As, Aui, Aid, Vu, Vi, Vd, A_layer

def get_f(t, pth='data/discharge_gb_box.nc', getF=-1, boxMethod=1, cal_clim=-1):
    # Get fresh water discharge
    import netCDF4 as nc

    print 'Getting freshwater discharge...'

    if getF==1:
        from gb_toolbox.gb_discharge import get_discharge_avgbox
        datapth = '/Users/chuning/projects/glacierbay/python/data/discharge_gb.nc'
        savepth = '/Users/chuning/Dropbox/projects/gb_boxmodel/data/' 
        get_discharge_avgbox(datapth, savepth, boxMethod=boxMethod)

    fh = nc.Dataset(pth, 'r')
    t_dis = fh.variables['t'][:]
    F = fh.variables['discharge'][:]
    fh.close()

    if cal_clim == 1:
        # Use climatology
        print 'Calculating discharge climatology...'
        pyt = mdates.num2date(t)
        yrday = np.array([pyt[i].timetuple().tm_yday for i in range(t.size)])

        pyt_dis = mdates.num2date(t_dis)
        yrday_dis = np.array([pyt_dis[i].timetuple().tm_yday for i in range(t_dis.size)])

        F_c = np.zeros((t.size, 4))
        for i in range(t.size):
            F_c[i, :] = np.nanmean(F[yrday_dis==yrday[i], :], axis=0)

        F = F_c
    else:
        F = interp1d(t_dis, F.T)(t).T
    return F

def get_wind(t, pth='data/juneau.csv', cal_clim=-1, filtwindow=15):
    # Get wind speed data from Juneau International Airport Station
    from box_gb.box_cdo import rd_cdo
    from gb_toolbox.buoy import smooth2a
    from scipy.interpolate import interp1d
    import matplotlib.dates as mdates

    print 'Getting wind speed from Juneau International Airport Station...'
    juneau = rd_cdo(pth)
    msk = juneau['awnd'].mask
    awnd = juneau['awnd'][~msk]
    # pdb.set_trace()
    pyt_awnd = juneau['pyt'][~msk]
    t_awnd = mdates.date2num(pyt_awnd)

    print 'Filtering dataset with a window of '+"%2d"%(2*filtwindow+1)+' days...'
    awnd = smooth2a(awnd, filtwindow)
    
    if cal_clim==1:
        print 'Calculating wind speed climatology...'
        pyt = mdates.num2date(t)
        yrday = np.array([pyt[i].timetuple().tm_yday for i in range(t.size)])
        yrday_awnd = np.array([pyt_awnd[i].timetuple().tm_yday for i in range(t_awnd.size)])

        awnd_c = np.zeros(t.size)
        for i in range(t.size):
            awnd_c[i] = np.nanmean(awnd[yrday_awnd==yrday[i]])
        awnd = awnd_c
    else:
        awnd = interp1d(t_awnd, awnd)(t)
    return awnd

def get_pacific(mt, var='s', lat_sp=57.875, lon_sp=-136.625, d_lim=[50, 100], uselowess=-1, filterwindow=60):
    # Get Pacific salinity
    import netCDF4 as nc

    # lat_sp = 57.875
    # lon_sp = -136.625
    # d_lim = [20, 100]

    # Load WOD salinity data
    fh = nc.Dataset('./data/ctd_wod.nc', 'r')
    mm = fh.variables['mm'][:]
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    d = fh.variables['d'][:]
    s = fh.variables[var][:]

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

    # Interpolate onto mt
    # Since num2date cannot deal with mt<1, add 366 days (1 year) if mt<1
    mt2 = np.zeros(mt.size)
    mt2[:] = mt
    mt2[mt2<1] = mt2[mt2<1]+366
    pyt = np.array(mdates.num2date(mt2))
    yd_mt = np.array([pyt[i].timetuple().tm_yday for i in range(pyt.size)])

    if uselowess == 1:
        # Use LOWESS filter
        import statsmodels.api as sm
        lowess = sm.nonparametric.lowess(sp_i, yd_i, frac=0.2)
        sp = interp1d(lowess[:, 0], lowess[:, 1])(yd_mt)

    else:
        # Use Butterworth and filtfilt

        from scipy.signal import butter, filtfilt

        fs = 1./(yd_i[1]-yd_i[0])  # Sampling frequency [day^-1]
        cutoff = 1./filterwindow  # Cutoff frequency [day^-1]
        b, a = butter(5, cutoff/(0.5*fs))

        # pdb.set_trace()

        sp_i = filtfilt(b, a, sp_i)
        sp = interp1d(yd_i, sp_i)(yd_mt)
    return sp

# -----------------------------------------------------------------------------------------
# working directory
rtdir = './'

# IC and BC data directory
icdir = './data/'

# Key parameters
getav = 1
getF = 0
cal_clim = 0
boxMethod = 1
F_scale = 1
tide_sn = 1
damp_mixing = 1
damp_wind_mixing = 1

# -----------------------------------------------------------------------------------------
# Constants
# Time step
t_init = mdates.date2num(datetime(1993, 01, 01))
t_end  = mdates.date2num(datetime(2008, 01, 01))
deltat = 1./2  # [day]

# Time
t = np.arange(t_init, t_end, deltat)

# Cg*beta*rho0*S0 [See Wang 2015 Msc Thesis]
Cbrs = 1e6  # [m^3s^-1]

# Relaxation Timescale [See Li et al 1999]
tr = 20.  # [day]

# Reference Salinity
S0 = 31.  # [PSU]

# Pacific Salinity
Spbar = 32.  # [PSU]

# Initial Salinity
Sinitbar = S0*np.ones(12)  # [PSU]

# Stratification Damp
deltaSbar_ui = np.array([2., 2., 1., 0.3])  # [PSU]
deltaSbar_id = np.array([0.1, 0.1, 0.1, 0.1])  # [PSU]

# Maximum Depth [m]
depth = 300

# Layer Thicknesses [m]
hu = 20
hi = 30
hd = depth-hu-hi

# Areas [m^2] and Volumes [m^3]
if getav == 1:
    # Calculate areas and volumes with digtal elevation model (DEM) outputs
    As, Aui, Aid, Vu, Vi, Vd = get_area_volume(hu, hi, hd, pth=icdir+'ARDEMv2.0.nc', boxMethod=boxMethod)[0:-1]
else:
    # Use a specified value
    # Here use cal_area_volume outputs to save computational time
    As  = np.array([3.12, 1.06, 3.56, 1.70])*1e8
    Aui = np.array([3.07, 0.96, 3.32, 1.37])*1e8
    Aid = np.array([2.65, 0.75, 2.74, 2.88])*1e8

    Vu = np.array([6.21,   2.06,  6.92,  3.13])*1e9
    Vi = np.array([14.38,  4.30, 15.40,  4.84])*1e9
    Vd = np.array([48.83,  7.66, 26.10,  0.50])*1e9

if boxMethod==2:
    A23u = 50.*1e3*hu
    A23i = 50.*1e3*hi

# Mixing Velocity
omega_ui = np.array([0.02, 0.02, 0.2, 50.])*1e-5  # [m s^-1]
omega_id = np.array([5., 5., 5., 3.])*1e-5  # [m s^-1]

if boxMethod==2:
    omega_23u = 1e-2  # [m s^-1]
    omega_23i = 1e-2  # [m s^-1]

# -----------------------------------------------------------------------------------------
# Freshwater Discharge
F = get_f(t, pth=icdir+'discharge_gb_box.nc', getF=getF, boxMethod=boxMethod, cal_clim=cal_clim)  # m^3day^-1
F = F/24/60/60  # m^3s^-1
F = F*F_scale

# Wind Speed (Juneau International Airport)
if damp_wind_mixing==1:
    awnd = get_wind(t, pth=icdir+'juneau.csv', cal_clim=cal_clim)
    awnd_crit = 6  # [m s^-1]
# -----------------------------------------------------------------------------------------
# Non-dimensionallize
Sp = Spbar/S0
deltaSui = deltaSbar_ui/S0
deltaSid = deltaSbar_id/S0

V = np.sum(Vu.sum()+Vi.sum()+Vd.sum())
ts = V/Cbrs/24/60/60
lambda_u = Vu/V
lambda_i = Vi/V
lambda_d = Vd/V

Rui = omega_ui*Aui/Cbrs
Rid = omega_id*Aid/Cbrs

if boxMethod==2:
    R23u = omega_23u*A23u/Cbrs
    R23i = omega_23i*A23i/Cbrs

if boxMethod==1:
    R12 = 0.1  # Cbrs23/Cbrs03
    R23 = 0.8  # Cbrs34/Cbrs03
    R3p = 0.5  # Cbrs4p/Cbrs03
elif boxMethod==2:
    R12 = 1  # Cbrs12/Cbrs03
    R2p = 1  # Cbrs2p/Cbrs03
    R3p = 1  # Cbrs3p/Cbrs03

Rt = 1./tr

f = F/Cbrs

