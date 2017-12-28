"""
Box model for Glacer Bay
The box model is developed based on the SoG box model.

Chuning Wang
"""

# ------------------------- import modules ------------------------------------
import ConfigParser
import os
import fnmatch
import glob
import csv
from datetime import datetime, timedelta

import numpy as np
from scipy import interpolate
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import path
from mpl_toolkits.basemap import Basemap
from scipy.signal import filtfilt


# ------------------------- class Box -----------------------------------------
class Box(object):
    """ class to store some basic info about the box model. """

    def __init__(self, info):
        self.info = info
        if 'box_method' not in self.info.keys():
            self.info['box_method'] = 1

        if 'sl_river' not in self.info.keys():
            self.info['sl_river'] = 'l'

        if 'sl_sp' not in self.info.keys():
            self.info['sl_sp'] = 'l'

        if 'clim' not in self.info.keys():
            self.info['clim'] = False

        if 'compare' not in self.info.keys():
            self.info['compare'] = False

        self.boxes = {}
        self.volumes = {}
        self.areas = {}
        self.consts = {}
        self.measurements = {}

        self.get_consts()

        if self.info['clim']:
            self.info['t0'] = 1.
            self.info['t1'] = 366.

        self.time = np.arange(self.info['t0'],
                              self.info['t1'],
                              self.consts['dt'])

        return None

    def __call__(self):
        self.get_box()
        self.get_area_volume()
        if self.info['box_method'] == 1:
            self._fix_box_method1()
        if self.info['box_method'] == 2:
            self._fix_box_method2()

        self.get_consts2()

        # get forcing data
        self.BoxRivers = BoxRivers(self)
        self.BoxRivers()
        self.BoxSp = BoxSp(self)
        self.BoxSp()

        if self.info['compare']:
            self.get_measurements()

        return None

    def get_box(self):
        """ generate a python dict contains the coordinates of the boxes. """

        if self.info['box_method'] == 1:
            print('Generating Box Type 1...')
            self.boxes['box0'] = np.array([
                [-136.65, 58.55],
                [-137.30, 58.65],
                [-137.30, 59.15],
                [-135.70, 59.15],
                [-135.60, 58.75],
                [-135.60, 58.55]])

            self.boxes['box1'] = np.array([
                [-136.65, 58.55],
                [-135.60, 58.55],
                [-135.60, 58.50],
                [-136.00, 58.35]])

            self.boxes['box2'] = np.array([
                [-136.65, 58.55],
                [-136.75, 58.25],
                [-136.50, 57.95],
                [-135.20, 57.95],
                [-135.40, 58.55],
                [-135.60, 58.55],
                [-135.60, 58.50],
                [-136.00, 58.35]])

        if self.info['box_method'] == 2:
            print('Generating Box Type 2...')
            self.boxes['box0'] = np.array([
                [-136.65, 58.65],
                [-137.30, 58.65],
                [-137.30, 59.15],
                [-136.40, 59.15],
                [-136.20, 58.75]])

            self.boxes['box1'] = np.array([
                [-136.20, 58.75],
                [-136.40, 59.15],
                [-135.70, 59.15],
                [-135.60, 58.75]])

            self.boxes['box2'] = np.array([
                [-136.65, 58.55],
                [-136.65, 58.65],
                [-136.20, 58.75],
                [-135.60, 58.75],
                [-135.60, 58.55]])

            self.boxes['box3'] = np.array([
                [-136.65, 58.55],
                [-135.60, 58.55],
                [-135.60, 58.50],
                [-136.00, 58.35]])

            self.boxes['box4'] = np.array([
                [-136.65, 58.55],
                [-136.75, 58.25],
                [-136.50, 57.95],
                [-135.20, 57.95],
                [-135.40, 58.55],
                [-135.60, 58.55],
                [-135.60, 58.50],
                [-136.00, 58.35]])

        self.boxes_list = self.boxes.keys()
        self.boxes_list.sort()

        return None

    def get_consts(self):
        """ parse model constants from box.in to self.consts. """

        config = ConfigParser.ConfigParser()
        config.read('box.in')
        for opt in config.options(self.info['case']):
            self.consts[opt.lower()] = float(config.get(self.info['case'], opt))

        return None

    def get_consts2(self):
        """ add some combination of constants """

        for boxi in self.boxes_list:
            idx = boxi[-1]
            try:
                self.consts['c' + idx + 'brs'] = \
                    self.consts['c' + idx] * self.consts['beta'] * \
                    self.consts['rho0'] * self.consts['s0'] * \
                    self.consts['s2d']
            except:
                pass

            try:
                self.consts['mui' + idx] = \
                    self.consts['omegaui' + idx] * self.areas['Aui' + idx]
                self.consts['mid' + idx] = \
                    self.consts['omegaid' + idx] * self.areas['Aid' + idx]
            except:
                pass

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
        depth_d = depth - self.consts['h_id']
        depth_d[depth_d < 0] = 0
        depth_i = depth - self.consts['h_si']
        depth_i[depth_i < 0] = 0

        volume_d = area*(depth_d)
        volume_i = area*(depth_i) - volume_d
        volume_s = volume - volume_i - volume_d

        area_id = area*(depth_d > 0)
        area_si = area*(depth_i > 0)
        area_s = area*(depth > 0)

        # calculate box total volume
        for box in self.boxes_list:
            idx = box[-1]
            box_i = box_idx == float(idx)
            self.volumes['V' + idx] = np.sum(volume*box_i)
            self.volumes['Vu' + idx] = np.sum(volume_s*box_i)
            self.volumes['Vi' + idx] = np.sum(volume_i*box_i)
            self.volumes['Vd' + idx] = np.sum(volume_d*box_i)
            self.areas['Au' + idx] = np.sum(area_s*box_i)
            self.areas['Aui' + idx] = np.sum(area_si*box_i)
            self.areas['Aid' + idx] = np.sum(area_id*box_i)

        return None

    def _fix_box_method1(self):
        """ combine intermediate and deep layer for box 1. """

        self.volumes['Vi1'] = self.volumes['Vi1'] + self.volumes['Vd1']
        self.volumes['Vd1'] = 0
        self.volumes['Vi2'] = self.volumes['Vi2'] + self.volumes['Vd2']
        self.volumes['Vd2'] = 0

        return None

    def _fix_box_method2(self):
        """ combine intermediate and deep layer for box 2. """

        self.volumes['Vd2'] = self.volumes['Vd2'] + self.volumes['Vd0'] + \
            self.volumes['Vd1']
        self.volumes['Vd0'] = 0
        self.volumes['Vd1'] = 0
        self.volumes['Vi3'] = self.volumes['Vi3'] + self.volumes['Vd3']
        self.volumes['Vd3'] = 0
        self.volumes['Vi4'] = self.volumes['Vi4'] + self.volumes['Vd4']
        self.volumes['Vd4'] = 0

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

        for box in self.boxes_list:
            lon = self.boxes[box][:, 0]
            lat = self.boxes[box][:, 1]
            xpos, ypos = m_h(lon, lat)
            m_h.plot(xpos, ypos, '--k')
            m_h.plot([xpos[0], xpos[-1]], [ypos[0], ypos[-1]], '--k')

        return None

    def get_measurements(self):
        """ get in-situ measurements climatology. """

        fid = open('./data/gb_salt_clim_' + str(self.info['box_method']) + '.txt')
        text = fid.readlines()
        fid.close()
        text = [x.strip().split(',') for x in text]
        for line in text:
            self.measurements[line[0]] = np.array(line[1:]).astype(float)

# ------------------------- class BoxReaders ----------------------------------


class BoxRivers(object):
    """ Class to parse river discharge data. """

    def __init__(self, Box):
        self.info = Box.info
        self.boxes = Box.boxes
        self.boxes_list = Box.boxes_list
        self.time = Box.time
        self.scale_factor = Box.consts['scale_river']
        if 'sl_rivers' not in self.info.keys():
            self.info['sl_rivers'] = 'l'

        self.data = {}
        if self.info['clim']:
            self.climatology = {}

    def __call__(self):
        if self.info['sl_rivers'] == 's':
            print('Loading freshwater discharge...')
            self.get_box_rivers()
        elif self.info['sl_rivers'] == 'l':
            print('Load river discharge data from netCDF file...')
            self.load_box_rivers()

        if self.scale_factor != 1:
            # scale river discharge
            self.scale_river()

        if self.info['clim']:
            self.cal_clim()
            self.data = self.climatology

        self.interp()

    def get_box_rivers(self):
        """ calculate river discharge in each box """

        # Load discharge data
        fin = nc.Dataset(self.info['river_raw_file_name'], 'r')
        self.data['time'] = fin.variables['t'][:]
        lat = fin.variables['lat'][:]
        lon = fin.variables['lon'][:]

        # # time conversion
        # don't do time conversion here - use ROMS convention 1900-01-01
        # dtime = (datetime(1900, 1, 1) - datetime(1, 1, 1)).days
        # self.data['time'] = self.data['time'] + dtime

        hydro_box = _get_box_idx(self.boxes, lon, lat)
        # Find coastal cells
        coast = fin.variables['coast'][:]
        hydro_box[coast.mask] = -1

        # Divide GB into several hydro regions
        print('Averaging data within each box...')

        for box in self.boxes_list:
            box_idx = box[-1]
            self.data['river' + box_idx] = np.zeros(self.data['time'].shape)

        for i in range(len(self.data['time'])):
            rslice = fin.variables['discharge'][i, :, :]
            rslice[rslice <= 0] = np.NaN
            for box in self.boxes_list:
                box_idx = box[-1]
                self.data['river' + box_idx][i] = np.nansum(
                    rslice[hydro_box == float(box_idx)])

        fin.close()

        # save data
        print('Saving data to ' + self.info['river_file_name'])
        fout = nc.Dataset(self.info['river_file_name'], 'w')
        fout.description = 'Glacier Bay freshwater discharge and deglaciation, sum of each box'

        fout.createDimension('time')

        t_nc = fout.createVariable('time', 'f8', ('time'))
        t_nc[:] = self.data['time']
        t_nc.units = 'days since 1900-01-01'

        for box in self.boxes_list:
            box_idx = box[-1]
            var_nc = fout.createVariable('river' + box_idx, 'f8', ('time'))
            var_nc[:] = self.data['river' + box_idx]
            var_nc.units = 'm3day-1'

        return None

    def load_box_rivers(self):
        ''' load data from netCDF file. '''

        fin = nc.Dataset(self.info['river_file_name'], 'r')
        self.data['time'] = fin.variables['time'][:]

        for box in self.boxes_list:
            box_idx = box[-1]
            self.data['river' + box_idx] = fin.variables['river' + box_idx][:]
        fin.close()

        return None

    def scale_river(self):
        """ scale the river discharge with a scaling factor defined in box.in. """

        for box in self.boxes_list:
            box_idx = box[-1]
            self.data['river' + box_idx] = \
                self.data['river' + box_idx]*self.scale_factor
        return None

    def cal_clim(self):
        """ calculate climatology. """

        dtime = nc.num2date(self.data['time'], 'days since 1900-01-01')
        self.data['yearday'] = np.array([i.timetuple().tm_yday for i in dtime])

        self.climatology['time'] = np.arange(366) + 1
        for box in self.boxes_list:
            box_idx = box[-1]
            self.climatology['river' + box_idx] = []

            for i in range(366):
                msk = self.data['yearday'] == i + 1
                self.climatology['river' + box_idx].append(self.data['river' + box_idx][msk].mean())
            self.climatology['river' + box_idx] = np.array(self.climatology['river' + box_idx])

        return None

    def interp(self):
        """ interpolate data onto given times. """

        for box in self.boxes_list:
            box_idx = box[-1]
            self.data['river' + box_idx] = np.interp(self.time,
                self.data['time'], self.data['river' + box_idx])

        self.data['time'] = self.time

        return None

    def plt_rivers(self):
        """ plot river discharge in each box. """

        plt.figure()
        for box in self.boxes_list:
            box_idx = box[-1]
            plt.plot(self.data['time'], self.data['river' + box_idx])

        plt.legend(self.boxes_list)

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
        print('Loading Juneau weather station data...')
        self.read_cdo()

    def read_cdo(self):
        """ read CDO csv dataset. """

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
    """
    class to parse Pacific salinity data.
    Only run this piece of code on Alnilam.

    2019/11/14 implement a new way to represent Sp.
    Use WOD measurements instead of SODA. SODA is biased
    high in my case.
    """

    def __init__(self, Box):
        self.info = Box.info
        self.time = Box.time
        self.data = {}
        if self.info['clim']:
            self.climatology = {}

        # if Box.info['compare']:
        #     Box.get_measurements()
        # self.measurements = Box.measurements

    def __call__(self):
        # if self.info['sl_sp'] == 's':
        #     print('Read Sp data from SODA files...')
        #     self.read_sp_soda()
        # else:
        #     print('Load Sp data from netCDF file...')
        #     self.load_sp()
        # self.filter()
        # self.interp_daily()
        # self.interp()

        self.get_sp_data()

        if self.info['clim']:
            self.cal_clim()
            # self.data = self.climatology

        self.interp()

    def get_sp_data(self):
        """ get sp data from measurements in Icy Strait. """

        try:
            fid = open('./data/gb_icy_strait.txt', 'r')
        except:
            fid = open('../data/gb_icy_strait.txt', 'r')
        text = fid.readlines()
        fid.close()
        text = [x.strip().split(',') for x in text]
        self.data['time'] = np.array(text[0][1:]).astype(float)
        self.data['salt'] = np.array(text[2][1:]).astype(float)
        self.data['temp'] = np.array(text[4][1:]).astype(float)

        # spline interpolate to daily
        time = np.arange(366) + 1
        spl = interpolate.splrep(self.data['time'], self.data['salt'])
        self.data['salt'] = interpolate.splev(time, spl)
        spl = interpolate.splrep(self.data['time'], self.data['temp'])
        self.data['temp'] = interpolate.splev(time, spl)
        self.data['time'] = time

        # interpolate to a longer time grid (30 years)
        time =  np.arange(366*10) + 1
        dtime = nc.num2date(time, 'days since 1900-01-01')
        self.data['yearday'] = np.array([i.timetuple().tm_yday for i in dtime])
        self.data['salt'] = np.interp(self.data['yearday'], self.data['time'], self.data['salt'])
        self.data['temp'] = np.interp(self.data['yearday'], self.data['time'], self.data['temp'])
        self.data['time'] = time

    def read_sp_soda(self, lon_lim=[-138.30+360, -136.30+360], lat_lim=[57.50, 58.50], z_lim=[150, 300]):
        """ read SODA data to get process the salinity profiles """

        # parse lon/lat range
        if 'sp_lon_lim' in self.info.keys():
            lon_lim = self.info['sp_lon_lim']
        else:
            self.info['sp_lon_lim'] = lon_lim
        if 'sp_lat_lim' in self.info.keys():
            lat_lim = self.info['sp_lat_lim']
        else:
            self.info['sp_lat_lim'] = lat_lim
        if 'sp_z_lim' in self.info.keys():
            z_lim = self.info['sp_z_lim']
        else:
            self.info['sp_z_lim'] = z_lim

        fgrid = _recursive_glob(self.info['sp_soda_dir_name'], pattern='*_grid.nc')
        flist = _recursive_glob(self.info['sp_soda_dir_name'], pattern='soda3.3.1_5dy_ocean_*.nc')
        self.data['time'] = []
        self.data['salt'] = []
        self.data['temp'] = []

        fin = nc.Dataset(fgrid[0])
        lon0 = fin.variables['geolon_t'][0, :]
        lat0 = fin.variables['geolat_t'][:, 0]
        z0 = fin.variables['st_ocean'][:]
        fin.close()
        msk_lon = (lon0 >= lon_lim[0]) & (lon0 <= lon_lim[1])
        msk_lat = (lat0 >= lat_lim[0]) & (lat0 <= lat_lim[1])
        msk_z = (z0 >= z_lim[0]) & (z0 <= z_lim[1])

        for fname in flist:
            fin = nc.Dataset(fname, 'r')
            self.data['time'].append(fin.variables['time'][0])
            tunit = fin.variables['time'].units
            salt_pr = fin.variables['salt'][0, msk_z, msk_lat, msk_lon]
            temp_pr = fin.variables['temp'][0, msk_z, msk_lat, msk_lon]
            # salt_pr = salt_pr.mean(axis=(1, 2))
            # temp_pr = temp_pr.mean(axis=(1, 2))
            salt_pr = salt_pr.mean(axis=1).mean(axis=1)
            temp_pr = temp_pr.mean(axis=1).mean(axis=1)
            self.data['salt'].append(salt_pr.mean())
            self.data['temp'].append(temp_pr.mean())

        # convert time from 1980 based to 1900 based
        self.data['time'] = nc.date2num(nc.num2date(self.data['time'] , tunit), 'days since 1900-01-01 00:00:00')

        # save to netCDF file
        fout = nc.Dataset(self.info['sp_file_name'], 'w')

        fout.createDimension('time')

        for var in ['time', 'salt', 'temp']:
            fout.createVariable(var, 'f8', ('time'))
            fout.variables[var][:] = self.data[var]

        fout.close()

        return None

    def load_sp(self):
        """ load data from netCDF file """

        fin = nc.Dataset(self.info['sp_file_name'])
        for var in ['time', 'salt', 'temp']:
            self.data[var] = fin.variables[var][:]
        fin.close()

        return None

    def interp_daily(self):
        """ interpolate data to daily. """

        time = np.arange(self.data['time'][0], self.data['time'][-1])
        self.data['salt'] = np.interp(time, self.data['time'], self.data['salt'])
        self.data['temp'] = np.interp(time, self.data['time'], self.data['temp'])
        self.data['time'] = time

        return None

    def cal_clim(self):
        """ calculate climatology. """

        dtime = nc.num2date(self.data['time'], 'days since 1900-01-01')
        self.data['yearday'] = np.array([i.timetuple().tm_yday for i in dtime])

        self.climatology['time'] = np.arange(366) + 1
        self.climatology['salt'] = []
        self.climatology['temp'] = []

        for i in range(366):
            msk = self.data['yearday'] == i + 1
            self.climatology['salt'].append(self.data['salt'][msk].mean())
            self.climatology['temp'].append(self.data['temp'][msk].mean())
        self.climatology['temp'] = np.array(self.climatology['temp'])
        self.climatology['salt'] = np.array(self.climatology['salt'])

        return None

    def interp(self):
        """ interpolate data onto given times. """

        self.data['salt'] = np.interp(self.time,
            self.data['time'], self.data['salt'])
        self.data['temp'] = np.interp(self.time,
            self.data['time'], self.data['temp'])

        self.data['time'] = self.time

        return None

    def filter(self, nfilter=5):
        """ filter Sp and Tp. """

        bfilt = np.ones(nfilter)/float(nfilter)
        afilt = 1
        self.data['salt'] = filtfilt(bfilt, afilt, self.data['salt'])
        self.data['temp'] = filtfilt(bfilt, afilt, self.data['temp'])

    # ------------------------- obsolete --------------------------------------------------
    # def _read_sp_wod(self, lon_lim=[-138.30, -136.30], lat_lim=[57.50, 58.50], depth=500):
    #     """ read WOD raw file and process the salinity profiles """

    #     from wodpy import wod

    #     # parse lon/lat range
    #     if 'sp_lon_lim' in self.info.keys():
    #         lon_lim = self.info['sp_lon_lim']
    #     else:
    #         self.info['sp_lon_lim'] = lon_lim
    #     if 'sp_lat_lim' in self.info.keys():
    #         lat_lim = self.info['sp_lat_lim']
    #     else:
    #         self.info['sp_lat_lim'] = lat_lim

    #     zlev = np.arange(depth)
    #     lat = []
    #     lon = []
    #     spr = []
    #     tpr = []
    #     time = []
    #     yearday = []

    #     fid = open(self.info['sp_raw_file_name'])
    #     cts = 0  # counter
    #     loop = True
    #     while loop:
    #         prof = wod.WodProfile(fid)
    #         lon0, lat0 = prof.longitude(), prof.latitude()

    #         if (lon0 > lon_lim[0]) & (lon0 < lon_lim[1]) & (lat0 > lat_lim[0]) & (lat0 < lat_lim[1]):
    #             zpr = prof.z()
    #             if len(zpr) > 1:
    #                 salt = prof.s()
    #                 if not np.all(salt.mask):
    #                     lat.append(lat0)
    #                     lon.append(lon0)
    #                     time.append(nc.date2num(prof.datetime(), 'days since 1900-01-01'))
    #                     yearday.append(prof.datetime().timetuple().tm_yday)

    #                     temp = prof.s()
    #                     zmsk2 = ~salt.mask
    #                     zpr = zpr[zmsk2]
    #                     salt = salt[zmsk2]
    #                     temp = temp[zmsk2]

    #                     zmin, zmax = zpr[0], zpr[-1]
    #                     zmsk = (zlev >= zmin) & (zlev <= zmax)
    #                     spr_i = np.zeros(depth)*np.NaN
    #                     spr_i[zmsk] = np.interp(zlev[zmsk], zpr, salt)
    #                     tpr_i = np.zeros(depth)*np.NaN
    #                     tpr_i[zmsk] = np.interp(zlev[zmsk], zpr, salt)
    #                     spr.append(spr_i)
    #                     tpr.append(tpr_i)

    #         loop = not prof.is_last_profile_in_file(fid)
    #         cts += 1

    #     self.data_raw['lat'] = np.array(lat)
    #     self.data_raw['lon'] = np.array(lon)
    #     self.data_raw['time'] = np.array(time)
    #     self.data_raw['yearday'] = np.array(yearday)
    #     self.data_raw['salt'] = np.array(spr)
    #     self.data_raw['temp'] = np.array(tpr)
    #     self.data_raw['depth'] = zlev

    #     # save to netCDF file
    #     fout = nc.Dataset(self.info['sp_file_name'], 'w')

    #     fout.createDimension('time')
    #     fout.createDimension('z', depth)

    #     for var in ['time', 'yearday', 'lon', 'lat']:
    #         fout.createVariable(var, 'f8', ('time'))
    #         fout.variables[var][:] = self.data_raw[var]

    #     for var in ['salt', 'temp']:
    #         fout.createVariable(var, 'f8', ('time', 'z'))
    #         fout.variables[var][:] = self.data_raw[var]

    #     fout.createVariable('depth', 'f8', ('z'))
    #     fout.variables['depth'][:] = self.data_raw['depth']

    #     fout.close()

    #     return None

    # def _cal_clim(self):
    #     """ calclate climatology for a single profile """

    #     salt_avg = np.nanmean(self.data_raw['salt'], axis=0)
    #     temp_avg = np.nanmean(self.data_raw['temp'], axis=0)

    #     # define deep water S and T values
    #     salt_d = 34.0
    #     temp_d = 4.3

    #     # fill up NaNs in data
    #     salt_avg[-1] = salt_d
    #     temp_avg[-1] = temp_d
    #     msk = ~np.isnan(temp_avg)
    #     temp_avg = np.interp(self.data_raw['depth'],
    #                          self.data_raw['depth'][msk], temp_avg[msk])
    #     salt_avg = np.interp(self.data_raw['depth'],
    #                          self.data_raw['depth'][msk], salt_avg[msk])

    #     # filter
    #     bfilt = np.ones(20)/20.
    #     afilt = 1
    #     temp_avg = filtfilt(bfilt, afilt, temp_avg)
    #     salt_avg = filtfilt(bfilt, afilt, salt_avg)

    #     # give salinity and temperature some seasonal variability
    #     self.climatology['time'] = np.array([0., 31., 59., 90., 120., 151.,
    #         181., 212., 243., 273., 304., 334.]) + 14
    #     salt_clim_s = np.array([31, 31, 30.5, 30.5, 30.5, 30.5,
    #         30, 29.5, 29, 30, 31, 31])
    #     temp_clim_s = np.array([8, 8, 8, 9, 10, 10, 10, 11, 11, 11, 10, 9])

    #     self.climatology['salt'] = np.zeros((len(self.climatology['time']),
    #                                          len(self.data_raw['depth'])))
    #     self.climatology['temp'] = np.zeros((len(self.climatology['time']),
    #                                          len(self.data_raw['depth'])))
    #     self.climatology['depth'] = self.data_raw['depth']

    #     for i in range(len(self.data_raw['depth'])):
    #         self.climatology['salt'][:, i] = (salt_avg[i] - salt_d)/(salt_avg[0] - salt_d)*(salt_clim_s - salt_d) + salt_d
    #         self.climatology['temp'][:, i] = (temp_avg[i] - temp_d)/(temp_avg[0] - temp_d)*(temp_clim_s - temp_d) + temp_d

    #     return None

# ------------------------- other functionals ---------------------------------------------


def _recursive_glob(in_dir, pattern='*'):
    ''' Search recursively for files matching a specified pattern. '''

    matches = []
    for root, dirnames, filenames in os.walk(in_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

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
