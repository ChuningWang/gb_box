"""
Box model for Glacer Bay
The box model is developed based on the SoG box model.

Chuning Wang
"""

# ------------------------- import modules ------------------------------------
import pdb
import glob
import csv
from datetime import datetime, timedelta

import numpy as np
import netCDF4 as nc
import matplotlib.dates as mdates
from matplotlib import path
from mpl_toolkits.basemap import Basemap
from scipy.signal import filtfilt

# ------------------------- constants -----------------------------------------
# working directory
ROOT_DIR = './'

# grid directory
GRID_FILE = '~/Documents/gb_roms/grd/GlacierBay_hr_grd.nc'

# FORCING data directory
FORCING_DIR = './data/'

# Constants
H_SI = 15
H_ID = 50

# Switches
GET_AV = 1
GET_F = 0
CAL_CLIM = 0
BOX_METHOD = 1
F_SCALE = 1
TIDE_SN = 1
DAMP_MIXING = 1
DAMP_WIND_MIXING = 1

# ------------------------- class Box -----------------------------------------


class Box(object):
    """ class to store some basic info about the box model. """

    def __init__(self, info, H_SI=10, H_ID=50):
        self.info = info
        self.h_si = H_SI
        self.h_id = H_ID
        if 'box_method' not in self.info.keys():
            self.info['box_method'] = 1
        elif 'sl_river' not in self.info.keys():
            self.info['sl_river'] = 'l'
        self.boxes = {}
        self.volumes = {}
        self.areas = {}
        return None

    def __call__(self):
        self.get_box()
        self.get_area_volume()
        if self.info['box_method'] == 1:
            self._fix_box_method1()

    def get_box(self):
        """ generate a python dict contains the coordinates of the boxes. """

        if self.info['box_method'] == 1:
            print('Generating Box Type 1...')
            self.boxes['box0'] = np.array([
                [-137.30, 59.15],
                [-136.50, 59.15],
                [-136.15, 58.75],
                [-135.60, 58.75],
                [-135.60, 58.55],
                [-136.65, 58.55],
                [-137.30, 58.65]])

            self.boxes['box1'] = np.array([
                [-136.50, 59.15],
                [-135.70, 59.15],
                [-135.60, 58.75],
                [-136.15, 58.75]])

            self.boxes['box2'] = np.array([
                [-136.65, 58.55],
                [-135.60, 58.55],
                [-135.60, 58.50],
                [-136.00, 58.35]])

            self.boxes['box3'] = np.array([
                [-136.65, 58.55],
                [-136.60, 58.35],
                [-136.50, 57.95],
                [-135.20, 57.95],
                [-135.30, 58.50],
                [-135.60, 58.55],
                [-135.60, 58.50],
                [-136.00, 58.35]])

        if self.info['box_method'] == 2:
            print('Generating Box Type 2...')
            self.boxes['box0'] = np.array([
                [-137.7, 59.1],
                [-137.7, 59.25],
                [-136.6, 59.25],
                [-136.15, 58.725],
                [-136.65, 58.725],
                [-137.15, 58.65]])

            self.boxes['box1'] = np.array([
                [-136.6, 59.25],
                [-136.3, 59.25],
                [-135.7, 59.0],
                [-135.7, 58.725],
                [-136., 58.825],
                [-136.15, 58.725],
                [-136.3, 58.9]])

            self.boxes['box2'] = np.array([
                [-136.15, 58.725],
                [-136., 58.825],
                [-135.7, 58.725],
                [-135.65, 58.55],
                [-136.65, 58.55],
                [-136.65, 58.725]])

            self.boxes['box3'] = np.array([
                [-136.65, 58.55],
                [-135.65, 58.55],
                [-135.65, 58.45],
                [-136.0, 58.375]])

            self.boxes['box4'] = np.array([
                [-136.65, 58.55],
                [-136.6, 58.35],
                [-136.5, 58.00],
                [-135.2, 58.00],
                [-135.3, 58.50],
                [-135.65, 58.45],
                [-136.0, 58.375]])

        elif self.info['box_method'] == 3:
            print ('Generating Box Type 3...')
            self.boxes['box0'] = np.array([[-137.7, 59.1],
                                           [-137.7, 59.25],
                                           [-136.6, 59.25],
                                           [-136.15, 58.725],
                                           [-136.65, 58.725],
                                           [-137.15, 58.65]])

            self.boxes['box1'] = np.array([[-136.6, 59.25],
                                           [-136.3, 59.25],
                                           [-135.7, 59.0],
                                           [-135.7, 58.725],
                                           [-136., 58.825],
                                           [-136.15, 58.725],
                                           [-136.3, 58.9]])

            self.boxes['box2'] = np.array([[-136.15, 58.725],
                                           [-136., 58.825],
                                           [-135.7, 58.725],
                                           [-135.6, 58.425],
                                           [-135.975, 58.375],
                                           [-136.0, 58.575]])

            self.boxes['box3'] = np.array([[-136.65, 58.725],
                                           [-136.15, 58.725],
                                           [-136.0, 58.575],
                                           [-135.975, 58.375],
                                           [-136.65, 58.575]])

        return None

    def get_area_volume(self):
        """ get box area and volume. """

        fin = nc.Dataset(self.info['grid_file_name'], 'r')
        lat = fin.variables['lat_rho'][:]
        lon = fin.variables['lon_rho'][:]
        depth = fin.variables['h'][:]
        msk = fin.variables['mask_rho'][:]

        depth[depth < 0] = 0

        # get indices for each box
        box_idx = _get_box_idx(self.boxes, lon, lat)
        # mask land
        box_idx[msk == 0] = -1

        # calculate area of each grid
        area = 1./(fin.variables['pm'][:]*fin.variables['pn'][:])
        fin.close()
        # calculate volume of each grid
        volume = area*depth

        # volume for each vertical layer
        depth_d = depth - self.h_id
        depth_d[depth_d < 0] = 0
        depth_i = depth - self.h_si
        depth_i[depth_i < 0] = 0

        volume_d = area*(depth_d)
        volume_i = area*(depth_i) - volume_d
        volume_s = volume - volume_i - volume_d

        area_id = area*(depth_d > 0)
        area_si = area*(depth_i > 0)
        area_s = area*(depth > 0)

        # calculate box total volume
        box_list = self.boxes.keys()
        for box in box_list:
            idx = box[-1]
            box_i = box_idx == float(idx)
            self.volumes['V' + idx] = np.sum(volume*box_i)
            self.volumes['Vs' + idx] = np.sum(volume_s*box_i)
            self.volumes['Vi' + idx] = np.sum(volume_i*box_i)
            self.volumes['Vd' + idx] = np.sum(volume_d*box_i)
            self.areas['As' + idx] = np.sum(area_s*box_i)
            self.areas['Asi' + idx] = np.sum(area_si*box_i)
            self.areas['Aid' + idx] = np.sum(area_id*box_i)

        return None

    def _fix_box_method1(self):
        """ combine intermediate and deep layer for box 2. """

        self.volumes['Vi2'] = self.volumes['Vi2'] + self.volumes['Vd2']
        self.volumes['Vd2'] = 0
        return None

    def plt_box(self):
        ''' plot boxes on a basemap. '''

        lat_min = 57.75
        lat_max = 59.25
        lat_0 = 0.5 * (lat_min + lat_max)

        lon_min = -137.5
        lon_max = -135.0
        lon_0 = 0.5 * (lon_min + lon_max)

        m_h = Basemap(projection='merc',
                      llcrnrlon=lon_min, llcrnrlat=lat_min,
                      urcrnrlon=lon_max, urcrnrlat=lat_max,
                      lat_0=lat_0, lon_0=lon_0,
                      resolution='i')

        # draw background
        m_h.fillcontinents(color='lightgrey', alpha=0.5)
        m_h.drawmeridians(np.arange(lon_min, lon_max, 0.5),
                          labels=[0, 0, 0, 1], fontsize=6, linewidth=.2)
        m_h.drawparallels(np.arange(lat_min, lat_max, 0.25),
                          labels=[1, 0, 0, 0], fontsize=6, linewidth=.2)

        box_list = self.boxes.keys()
        for box in box_list:
            lon = self.boxes[box][:, 0]
            lat = self.boxes[box][:, 1]
            xpos, ypos = m_h(lon, lat)
            m_h.plot(xpos, ypos, '--k')
            m_h.plot([xpos[0], xpos[-1]], [ypos[0], ypos[-1]], '--k')

        return None

# ------------------------- class BoxReaders ----------------------------------


class BoxRivers(object):
    """ Class to parse river discharge data. """

    def __init__(self, Box):
        self.info = Box.info
        self.boxes = Box.boxes
        if 'sl_rivers' not in self.info.keys():
            self.info['sl_rivers'] = 'l'

        self.data = {}

    def __call__(self):
        if self.info['sl_rivers'] == 's':
            self.get_box_rivers()
        elif self.info['sl_rivers'] == 'l':
            self.load_box_rivers()

    def get_box_rivers(self):
        """ calculate river discharge in each box """

        # Load discharge data
        print('Loading freshwater discharge...')
        fin = nc.Dataset(self.info['river_raw_file_name'], 'r')
        self.data['time'] = fin.variables['t'][:]
        lat = fin.variables['lat'][:]
        lon = fin.variables['lon'][:]

        # time conversion
        dtime = (datetime(1900, 1, 1) - datetime(1, 1, 1)).days
        self.data['time'] = self.data['time'] + dtime

        hydro_box = _get_box_idx(self.boxes, lon, lat)
        # Find coastal cells
        coast = fin.variables['coast'][:]
        hydro_box[coast.mask] = -1

        fin.close()

        # Divide GB into several hydro regions
        print('Averaging data within each box...')

        box_list = self.boxes.keys()
        for box in box_list:
            box_idx = box[-1]
            self.data['river' + box_idx] = np.zeros(self.data['time'].shape)

        for i in range(len(self.data['time'])):
            rslice = fin.variables['discharge'][i, :, :]
            rslice[rslice <= 0] = np.NaN
            for box in box_list:
                box_idx = box[-1]
                self.data['river' + box_idx][i] = np.nansum(
                    rslice[hydro_box == float(box_idx)])

        # save data
        print 'Saving data to ' + self.info['river_file_name']
        fout = nc.Dataset(self.info['river_file_name'], 'w')
        fout.description = 'Glacier Bay freshwater discharge and deglaciation, sum of each box'

        fout.createDimension('time')

        t_nc = fout.createVariable('time', 'f8', ('time'))
        t_nc[:] = self.data['time']
        t_nc.units = 'days since 0001-01-01'

        for box in box_list:
            box_idx = box[-1]
            var_nc = fout.createVariable('river' + box_idx, 'f8', ('time'))
            var_nc[:] = self.data['river' + box_idx]
            var_nc.units = 'm3day-1'

        return None

    def load_box_rivers(self):
        ''' load data from netCDF file. '''

        print('Load data from netCDF file...')
        fin = nc.Dataset(self.info['river_file_name'], 'r')
        self.data['time'] = fin.variables['time'][:]

        box_list = self.boxes.keys()
        for box in box_list:
            box_idx = box[-1]
            self.data['river' + box_idx] = fin.variables['river' + box_idx][:]
        fin.close()

        return None


class BoxCDO(object):
    """
    class to parse Juneau wind data.

    Here is the documentation from CDO - 

    The five core values are:
    PRCP = Precipitation (mm or inches as per user preference, inches to hundredths on Daily Form pdf file)
    SNOW = Snowfall (mm or inches as per user preference, inches to tenths on Daily Form pdf file)
    SNWD = Snow depth (mm or inches as per user preference, inches on Daily Form pdf file)
    TMAX = Maximum temperature (Fahrenheit or Celsius as per user preference, Fahrenheit to tenths on Daily Form pdf file
    TMIN = Minimum temperature (Fahrenheit or Celsius as per user preference, Fahrenheit to tenths on Daily Form pdf file

    The other values are:
    ACMC = Average cloudiness midnight to midnight from 30-second ceilometer data (percent)
    ACMH = Average cloudiness midnight to midnight from manual observations (percent)
    ACSC = Average cloudiness sunrise to sunset from 30-second ceilometer data (percent)
    ACSH = Average cloudiness sunrise to sunset from manual observations (percent)
    AWND = Average daily wind speed (meters per second or miles per hour as per user preference)
    DAEV = Number of days included in the multiday evaporation total (MDEV)
    DAPR = Number of days included in the multiday precipitation total (MDPR)
    DASF = Number of days included in the multiday snowfall total (MDSF)
    DATN = Number of days included in the multiday minimum temperature (MDTN)
    DATX = Number of days included in the multiday maximum temperature (MDTX)
    DAWM = Number of days included in the multiday wind movement (MDWM)
    DWPR = Number of days with non-zero precipitation included in multiday precipitation total (MDPR)
    EVAP = Evaporation of water from evaporation pan (mm or inches as per user preference, or hundredths of inches on Daily Form pdf file)
    FMTM = Time of fastest mile or fastest 1-minute wind (hours and minutes, i.e., HHMM)
    FRGB = Base of frozen ground layer (cm or inches as per user preference)
    FRGT = Top of frozen ground layer (cm or inches as per user preference)
    FRTH = Thickness of frozen ground layer (cm or inches as per user preference)
    GAHT = Difference between river and gauge height (cm or inches as per user preference)
    MDEV = Multiday evaporation total (mm or inches as per user preference; use with DAEV)
    MDPR = Multiday precipitation total (mm or inches as per user preference; use with DAPR and DWPR, if available)
    MDSF = Multiday snowfall total (mm or inches as per user preference)
    MDTN = Multiday minimum temperature (Fahrenheit or Celsius as per user preference ; use with DATN)
    MDTX = Multiday maximum temperature (Fahrenheit or Celsius as per user preference ; use with DATX)
    MDWM = Multiday wind movement (miles or km as per user preference)
    MNPN = Daily minimum temperature of water in an evaporation pan (Fahrenheit or Celsius as per user preference)
    MXPN = Daily maximum temperature of water in an evaporation pan (Fahrenheit or Celsius as per user preference)

    PGTM = Peak gust time (hours and minutes, i.e., HHMM)
    PSUN = Daily percent of possible sunshine (percent)
    SN*# = Minimum soil temperature where * corresponds to a code
    for ground cover and # corresponds to a code for soil depth (Fahrenheit or Celsius as per user preference)
    Ground cover codes include the following: 0 = unknown
    1 = grass
    2 = fallow
    3 = bare ground 4 = brome grass 5 = sod
    6 = straw mulch 7 = grass muck 8 = bare muck
    Depth codes include the following: 1 = 5 cm
    2 = 10 cm 3 = 20 cm 4 = 50 cm 5 = 100 cm 6 = 150 cm 7 = 180 cm
    SX*# = Maximum soil temperature where * corresponds to a code for ground
    cover and # corresponds to a code for soil depth. See SN*# for depth codes. (Fahrenheit or
    Celsius as per user preference)
    THIC = Thickness of ice on water (inches or mm as per user preference)
    TOBS = Temperature at the time of observation (Fahrenheit or Celsius as per user preference)
    TSUN = Daily total sunshine (minutes)
    WDF1 = Direction of fastest 1-minute wind (degrees)
    WDF2 = Direction of fastest 2-minute wind (degrees)
    WDF5 = Direction of fastest 5-second wind (degrees)
    WDFG = Direction of peak wind gust (degrees)
    WDFI = Direction of highest instantaneous wind (degrees)
    WDFM = Fastest mile wind direction (degrees)
    WDMV = 24-hour wind movement (km or miles as per user preference, miles on Daily Form pdf file)
    WESD = Water equivalent of snow on the ground (inches or mm as per user preference)
    WESF = Water equivalent of snowfall (inches or mm as per user preference)
    WSF1 = Fastest 1-minute wind speed (miles per hour or meters per second as per user preference)
    WSF2 = Fastest 2-minute wind speed (miles per hour or meters per second as per user preference)
    WSF5 = Fastest 5-second wind speed (miles per hour or meters per second as per user preference)
    WSFG = Peak guest wind speed (miles per hour or meters per second as per user preference)
    WSFI = Highest instantaneous wind speed (miles per hour or meters per second as per user preference)
    WSFM = Fastest mile wind speed (miles per hour or meters per second as per user preference)
    WT** = Weather Type where ** has one of the following values:
    01 = Fog, ice fog, or freezing fog (may include heavy fog)
    02 = Heavy fog or heaving freezing fog (not always distinguished from fog)
    03 = Thunder
    04 = Ice pellets, sleet, snow pellets, or small hail
    05 = Hail (may include small hail)
    06 = Glaze or rime
    07 = Dust, volcanic ash, blowing dust, blowing sand, or blowing obstruction 08 = Smoke or haze
    09 = Blowing or drifting snow
    10 = Tornado, waterspout, or funnel cloud
    11 = High or damaging winds
    12 = Blowing spray
    13 = Mist
    14 = Drizzle
    15 = Freezing drizzle
    16 = Rain (may include freezing rain, drizzle, and freezing drizzle) 17 = Freezing rain
    18 = Snow, snow pellets, snow grains, or ice crystals
    19 = Unknown source of precipitation
    21 = Ground fog
    22 = Ice fog or freezing fog
    """

    def __init__(self, Box):
        self.info = Box.info
        self.data = {}

    def __call__(self):
        self.read_cdo()

    def read_cdo(self):
        """ read CDO csv dataset. """

        print 'Loading Juneau weather station data...'
        data = []
        header = []
        fin = open(self.info['wind_file_name'])
        reader = csv.reader(fin)
        for row in reader:
            if row[0][:7] == 'STATION':
                header.append(row)    # read header
                continue
            data.append(row)
        fin.close()

        data = np.asarray(data)
        header = header[0]

        for i, hdr in enumerate(header):
            if hdr == 'STATION':
                self.info['station'] = data[0, i]
            elif hdr == 'STATION_NAME':
                self.info['station_long_name'] = data[0, i]
            elif hdr == 'ELEVATION':
                self.data['elev'] = data[:, i].astype(np.float)
            elif hdr == 'LATITUDE':
                self.data['lat'] = data[:, i].astype(np.float)
            elif hdr == 'LONGITUDE':
                self.data['lon'] = data[:, i].astype(np.float)
            elif hdr == 'DATE':
                time = np.array([datetime.strptime(t_i, "%Y%m%d") for t_i in data[:, i]])
                self.data['time'] = nc.date2num(time, 'days since 1900-01-01')
            else:
                dslice = data[:, i].astype(np.float)
                dslice = np.ma.masked_where(dslice == -9999, dslice)
                self.data[hdr.lower()] = dslice

        return None


class BoxSp(object):
    """ class to parse Pacific salinity data. """

    def __init__(self, Box):
        self.info = Box.info
        if 'sl_sp' not in self.info.keys():
            self.info['sl_sp'] = 'l'
        self.data_raw = {}
        self.data = {}
        self.climatology = {}

    def __call__(self):
        if self.info['sl_sp'] == 's':
            self.read_sp()
        else:
            self.load_sp()

        self.cal_clim()
        self.cal_sp()

    def read_sp_wod(self, lon_lim=[-138.30, -136.30], lat_lim=[57.50, 58.50], depth=500):
        """ read WOD raw file and process the salinity profiles """

        from wodpy import wod

        # parse lon/lat range
        if 'sp_lon_lim' in self.info.keys():
            lon_lim = self.info['sp_lon_lim']
        else:
            self.info['sp_lon_lim'] = lon_lim
        if 'sp_lat_lim' in self.info.keys():
            lat_lim = self.info['sp_lat_lim']
        else:
            self.info['sp_lat_lim'] = lat_lim

        zlev = np.arange(depth)
        lat = []
        lon = []
        spr = []
        tpr = []
        time = []
        yearday = []

        fid = open(self.info['sp_raw_file_name'])
        cts = 0  # counter
        loop = True
        while loop:
            prof = wod.WodProfile(fid)
            lon0, lat0 = prof.longitude(), prof.latitude()

            if (lon0 > lon_lim[0]) & (lon0 < lon_lim[1]) & (lat0 > lat_lim[0]) & (lat0 < lat_lim[1]):
                zpr = prof.z()
                if len(zpr) > 1:
                    salt = prof.s()
                    if not np.all(salt.mask):
                        lat.append(lat0)
                        lon.append(lon0)
                        time.append(nc.date2num(prof.datetime(), 'days since 1900-01-01'))
                        yearday.append(prof.datetime().timetuple().tm_yday)

                        temp = prof.s()
                        zmsk2 = ~salt.mask
                        zpr = zpr[zmsk2]
                        salt = salt[zmsk2]
                        temp = temp[zmsk2]

                        zmin, zmax = zpr[0], zpr[-1]
                        zmsk = (zlev >= zmin) & (zlev <= zmax)
                        spr_i = np.zeros(depth)*np.NaN
                        spr_i[zmsk] = np.interp(zlev[zmsk], zpr, salt)
                        tpr_i = np.zeros(depth)*np.NaN
                        tpr_i[zmsk] = np.interp(zlev[zmsk], zpr, salt)
                        spr.append(spr_i)
                        tpr.append(tpr_i)

            loop = not prof.is_last_profile_in_file(fid)
            cts += 1

        self.data_raw['lat'] = np.array(lat)
        self.data_raw['lon'] = np.array(lon)
        self.data_raw['time'] = np.array(time)
        self.data_raw['yearday'] = np.array(yearday)
        self.data_raw['salt'] = np.array(spr)
        self.data_raw['temp'] = np.array(tpr)
        self.data_raw['depth'] = zlev

        # save to netCDF file
        fout = nc.Dataset(self.info['sp_file_name'], 'w')

        fout.createDimension('time')
        fout.createDimension('z', depth)

        for var in ['time', 'yearday', 'lon', 'lat']:
            fout.createVariable(var, 'f8', ('time'))
            fout.variables[var][:] = self.data_raw[var]

        for var in ['salt', 'temp']:
            fout.createVariable(var, 'f8', ('time', 'z'))
            fout.variables[var][:] = self.data_raw[var]

        fout.createVariable('depth', 'f8', ('z'))
        fout.variables['depth'][:] = self.data_raw['depth']

        fout.close()

        return None

    def read_sp_soda(self, lon_lim=[-138.30, -136.30], lat_lim=[57.50, 58.50], depth=500):
        """ read SODA data to get process the salinity profiles """

        from wodpy import wod

        # parse lon/lat range
        if 'sp_lon_lim' in self.info.keys():
            lon_lim = self.info['sp_lon_lim']
        else:
            self.info['sp_lon_lim'] = lon_lim
        if 'sp_lat_lim' in self.info.keys():
            lat_lim = self.info['sp_lat_lim']
        else:
            self.info['sp_lat_lim'] = lat_lim

        zlev = np.arange(depth)
        lat = []
        lon = []
        spr = []
        tpr = []
        time = []
        yearday = []

        fid = open(self.info['sp_raw_file_name'])
        cts = 0  # counter
        loop = True
        while loop:
            prof = wod.WodProfile(fid)
            lon0, lat0 = prof.longitude(), prof.latitude()

            if (lon0 > lon_lim[0]) & (lon0 < lon_lim[1]) & (lat0 > lat_lim[0]) & (lat0 < lat_lim[1]):
                zpr = prof.z()
                if len(zpr) > 1:
                    salt = prof.s()
                    if not np.all(salt.mask):
                        lat.append(lat0)
                        lon.append(lon0)
                        time.append(nc.date2num(prof.datetime(), 'days since 1900-01-01'))
                        yearday.append(prof.datetime().timetuple().tm_yday)

                        temp = prof.s()
                        zmsk2 = ~salt.mask
                        zpr = zpr[zmsk2]
                        salt = salt[zmsk2]
                        temp = temp[zmsk2]

                        zmin, zmax = zpr[0], zpr[-1]
                        zmsk = (zlev >= zmin) & (zlev <= zmax)
                        spr_i = np.zeros(depth)*np.NaN
                        spr_i[zmsk] = np.interp(zlev[zmsk], zpr, salt)
                        tpr_i = np.zeros(depth)*np.NaN
                        tpr_i[zmsk] = np.interp(zlev[zmsk], zpr, salt)
                        spr.append(spr_i)
                        tpr.append(tpr_i)

            loop = not prof.is_last_profile_in_file(fid)
            cts += 1

        self.data_raw['lat'] = np.array(lat)
        self.data_raw['lon'] = np.array(lon)
        self.data_raw['time'] = np.array(time)
        self.data_raw['yearday'] = np.array(yearday)
        self.data_raw['salt'] = np.array(spr)
        self.data_raw['temp'] = np.array(tpr)
        self.data_raw['depth'] = zlev

        # save to netCDF file
        fout = nc.Dataset(self.info['sp_file_name'], 'w')

        fout.createDimension('time')
        fout.createDimension('z', depth)

        for var in ['time', 'yearday', 'lon', 'lat']:
            fout.createVariable(var, 'f8', ('time'))
            fout.variables[var][:] = self.data_raw[var]

        for var in ['salt', 'temp']:
            fout.createVariable(var, 'f8', ('time', 'z'))
            fout.variables[var][:] = self.data_raw[var]

        fout.createVariable('depth', 'f8', ('z'))
        fout.variables['depth'][:] = self.data_raw['depth']

        fout.close()

        return None

    def load_sp(self):
        """ load data from netCDF file """

        fin = nc.Dataset(self.info['sp_file_name'])
        for var in ['time', 'yearday', 'lat', 'lon', 'depth', 'salt', 'temp']:
            self.data_raw[var] = fin.variables[var][:]
        fin.close()

        return None

    def cal_clim(self):
        """ calclate climatology for a single profile """

        salt_avg = np.nanmean(self.data_raw['salt'], axis=0)
        temp_avg = np.nanmean(self.data_raw['temp'], axis=0)

        # define deep water S and T values
        salt_d = 34.0
        temp_d = 4.3

        # fill up NaNs in data
        salt_avg[-1] = salt_d
        temp_avg[-1] = temp_d
        msk = ~np.isnan(temp_avg)
        temp_avg = np.interp(self.data_raw['depth'],
                             self.data_raw['depth'][msk], temp_avg[msk])
        salt_avg = np.interp(self.data_raw['depth'],
                             self.data_raw['depth'][msk], salt_avg[msk])

        # filter
        bfilt = np.ones(20)/20.
        afilt = 1
        temp_avg = filtfilt(bfilt, afilt, temp_avg)
        salt_avg = filtfilt(bfilt, afilt, salt_avg)

        # give salinity and temperature some seasonal variability
        self.climatology['time'] = np.array([0., 31., 59., 90., 120., 151.,
            181., 212., 243., 273., 304., 334.]) + 14
        salt_clim_s = np.array([31, 31, 30.5, 30.5, 30.5, 30.5,
            30, 29.5, 29, 30, 31, 31])
        temp_clim_s = np.array([8, 8, 8, 9, 10, 10, 10, 11, 11, 11, 10, 9])

        self.climatology['salt'] = np.zeros((len(self.climatology['time']),
                                             len(self.data_raw['depth'])))
        self.climatology['temp'] = np.zeros((len(self.climatology['time']),
                                             len(self.data_raw['depth'])))
        self.climatology['depth'] = self.data_raw['depth']

        for i in range(len(self.data_raw['depth'])):
            self.climatology['salt'][:, i] = (salt_avg[i] - salt_d)/(salt_avg[0] - salt_d)*(salt_clim_s - salt_d) + salt_d
            self.climatology['temp'][:, i] = (temp_avg[i] - temp_d)/(temp_avg[0] - temp_d)*(temp_clim_s - temp_d) + temp_d

        return None

    def cal_sp(self, zmin=150, zmax=300):
        """ calcualte sp from climatology data. """

        self.data['time'] = self.data_raw['time']
        self.data['yearday'] = self.data_raw['yearday']
        mskz = (self.data_raw['depth'] >= zmin) & (self.data_raw['depth'] <= zmax)
        self.data['sp'] = np.nanmean(self.data_raw['salt'][:, mskz], axis=1)

        return None

        # import statsmodels.api as sm

        # # Load WOD salinity data
        # fh = nc.Dataset('./data/s_wod.nc', 'r')
        # mm = fh.variables['mm'][:]
        # lat = fh.variables['lat'][:]
        # lon = fh.variables['lon'][:]
        # d = fh.variables['d'][:]
        # s = fh.variables['s'][:]

        # # Select data at coordinate (lon_sp, lat_sp)
        # msk1 = lat==lat_sp
        # msk2 = lon==lon_sp
        # msk3 = (d>=d_lim[0]) & (d<=d_lim[1])
        # sp = s[:, msk3, msk1, msk2]
        # sp = np.mean(sp, axis=1)

        # yd = np.zeros(12)
        # for i in range(yd.size):
        #     yd[i] = datetime(1, i+1, 15).timetuple().tm_yday

        # yd = np.concatenate((yd-366, yd, yd+366))
        # sp = np.concatenate((sp, sp, sp))

        # yd_i = np.arange(-300, 300*2)
        # sp_i = interp1d(yd,sp)(yd_i)

        # # LOWESS filter
        # lowess = sm.nonparametric.lowess(sp_i, yd_i, frac=0.2)
        # # Interpolate onto mt
        # # Since num2date cannot deal with mt<1, add 366 days (1 year) if mt<1
        # mt2 = np.zeros(mt.size)
        # mt2[:] = mt
        # mt2[mt2<1] = mt2[mt2<1]+366
        # pyt = np.array(mdates.num2date(mt2))
        # yd_mt = np.array([pyt[i].timetuple().tm_yday for i in range(pyt.size)])
        # sp = interp1d(lowess[:, 0], lowess[:, 1])(yd_mt)

        # # pdb.set_trace()

        # return sp

# ------------------------- other functionals ---------------------------------------------


def _get_rivers(data_path, save_path, lat_lim, lon_lim):
    """ load river discharge from raw discharge file.
    :param data_path: raw discharge data path
    :param save_path: output data path
    :param lat_lim: latitude range
    :param lon_lim: longitude range
    :return: None
    """

    # Reading lat lon
    print('Loading Lat & lon...')
    fin = nc.Dataset(data_path + 'lat_lon.nc')
    lat = fin.variables['lat'][:]
    lon = fin.variables['lon'][:]
    lon = lon - 360
    fin.close()

    # Cut out useful portion
    msk_gb = (lat > lat_lim[0]) & (lat < lat_lim[1]) & \
             (lon > lon_lim[0]) & (lon < lon_lim[1])
    msk_lat = np.any(msk_gb, axis=1)
    msk_lon = np.any(msk_gb, axis=0)

    lat = lat[msk_lat, :][:, msk_lon]
    lon = lon[msk_lat, :][:, msk_lon]

    # Load data
    print('Getting file names...')
    flist = glob.glob(data_path + 'discharge' + '*.nc')

    # Initiate
    print 'Loading ' + flist[0] + '...'
    fin = nc.Dataset(flist[0], 'r')
    dis = fin.variables['discharge'][:, :, msk_lat, :][:, :, :, msk_lon]
    time = fin.variables['time'][:]
    t_ini = datetime.strptime(fin.variables['time'].units[12:],
                              '%Y-%m-%d %H:%M:%S')
    pyt = np.array([t_ini + timedelta(hours=time[i]) for i in range(time.size)])

    for file_name in flist[1:]:
        print('Loading ' + file_name + '...')
        fh = nc.Dataset(file_name, 'r')
        d_in = fh.variables['discharge'][:]
        d_in = d_in[:, :, msk_lat, :][:, :, :, msk_lon]
        time = fh.variables['time'][:]
        t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
        pyt_in = np.array([t_ini + timedelta(hours=time[i]) for i in range(t.size)])
        fh.close()

        dis = np.concatenate([dis, d_in], axis=0)
        pyt = np.concatenate([pyt, pyt_in])

    # mask out invalid data
    print('Setting invalid data to NaN...')
    dis = np.squeeze(dis)
    dis[dis < -1000] = np.nan

    lon_dim = dis.shape[1]
    lat_dim = dis.shape[2]

    mtime = mdates.date2num(pyt)

    # Load coast cell file
    print('Getting coast cells...')
    fin = nc.Dataset(data_path + 'lat_lon.nc', 'r')
    coast_cells = np.squeeze(fin.variables['coast_cells'][:, msk_lat, :][:, :, msk_lon])
    fin.close()

    # write data into netCDF file
    print('Saving data as netCDF4 file...')
    fout = nc.Dataset(save_path + 'discharge_gb.nc', 'w')
    fout.description = 'Glacier Bay river discharge and deglaciation'

    fout.createDimension('time', None)
    fout.createDimension('lat', lat_dim)
    fout.createDimension('lon', lon_dim)

    t_nc = fout.createVariable('t', 'f8', ('time'))
    lat_nc = fout.createVariable('lat', 'f8', ('lon', 'lat'))
    lon_nc = fout.createVariable('lon', 'f8', ('lon', 'lat'))
    coast_nc = fout.createVariable('coast_cells', 'f8', ('lon', 'lat'))
    d_nc = fout.createVariable('discharge', 'f8', ('time', 'lon', 'lat'))

    t_nc[:] = mtime
    lat_nc[:, :] = lat
    lon_nc[:, :] = lon
    coast_nc[:, :] = coast_cells
    d_nc[:, :, :] = dis

    t_nc.units = 'days since 0001-01-01'
    d_nc.units = 'm^3day^-1'

    fout.close()

    return None


def _get_box_idx(boxes, lon, lat):
    """
    identify with box each grid point belongs to.
    : param boxes: boxes dict from get_box
    : param lon, lat: geo coordinates
    : return: box_idx
    """

    eta, xi = lon.shape
    box_idx = -1*np.ones(eta*xi)
    box_list = boxes.keys()
    path_points = {}
    for box in box_list:
        path_points[box] = path.Path(boxes[box])

    coords = []
    for i in range(eta):
        for j in range(xi):
            coords.append((lon[i, j], lat[i, j]))

    for box in box_list:
        idx = int(box[-1])
        box_i = path_points[box].contains_points(coords)
        box_idx[box_i] = idx

    box_idx = box_idx.reshape((eta, xi))

    return box_idx
