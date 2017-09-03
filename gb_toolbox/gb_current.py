import numpy as np
import urllib2
import xmltodict
import datetime as dt
from scipy.interpolate import interp1d
from scipy.signal import buttord, butter, filtfilt
import netCDF4 as nc

class get_noaa_current():
    def __init__(self, info):
        self.info = info
        if 'prod' not in self.info.keys():
            self.info['prod'] = 'currents'
        if 'units' not in self.info.keys():
            self.info['units'] = 'metric'
        if 'tz' not in self.info.keys():
            self.info['tz'] = 'gmt'
        if 'fmt' not in self.info.keys():
            self.info['fmt'] = 'csv'
        if 't0' not in self.info.keys():
            self.info['t0'] = '1900-01-01 00:00:00'
        if 'sl' not in self.info.keys():
            self.info['sl'] = 'l'
        if 'Wp_hrs' not in self.info.keys():
            self.info['Wp_hrs'] = -1

    def __call__(self):
        print 'Formating download urls'
        self.get_metadata()
        self.make_url()
        if self.info['sl']=='s':
            self.get_current()
        if self.info['sl']=='l':
            self.load_data()
        self.compute_uv()
        if self.info['Wp_hrs']>0:
            self.filter()
        if self.info['sl']=='s':
            self.save_data()

        return None

    def get_metadata(self):
        self.info['stn_url'] = 'http://tidesandcurrents.noaa.gov/mdapi/v0.6/webapi/stations/' + \
                               self.info['stn'] + '.xml'

        page = urllib2.urlopen(self.info['stn_url'])
        s_info = xmltodict.parse(page.read())
        s_info = s_info['Stations']['Station']
        self.info['stn_name'] = s_info['name']
        self.info['lat'] = float(s_info['lat'])
        self.info['lon'] = float(s_info['lng'])
        self.info['deployed'] = s_info['deployed']
        self.info['retrieved'] = s_info['retrieved']
        self.info['tz_offset'] = float(s_info['timezone_offset'])
        self.info['dply_url'] = s_info['deployments']['@self']
        self.info['bins_url'] = s_info['bins']['@self']

        page = urllib2.urlopen(self.info['dply_url'])
        d_info = xmltodict.parse(page.read())
        d_info = d_info['Deployments']
        self.info['depth'] = float(d_info['depth'])

        page = urllib2.urlopen(self.info['bins_url'])
        b_info = xmltodict.parse(page.read())
        b_info = b_info['Bins']['Bin']
        self.info['zlevs'] = len(b_info)
        z = np.array([float(b_info[i]['depth']) for i in range(self.info['zlevs'])])
        self.z = z

    def make_url(self):
        self.url = ['https://tidesandcurrents.noaa.gov/api/datagetter?' + \
                    'begin_date=' + self.info['bdate'] + \
                    '&end_date=' + self.info['edate'] + \
                    '&station=' + self.info['stn'] + \
                    '&product=' + self.info['prod'] + \
                    '&bin=' + str(bin) + \
                    '&units=' + self.info['units'] + \
                    '&time_zone=' + self.info['tz'] + \
                    '&application=web_services' + \
                    '&format=' + self.info['fmt'] for bin in range(1, self.info['zlevs']+1)]

        return None

    def wget_current(self, idx):
        page = urllib2.urlopen(self.url[idx])
        web = page.read().split('\n')
        web.pop(0)
        web.pop(-1)
        ct = len(web)
        for i in range(ct):
            web[i] = web[i].split(',')

        time = np.array([dt.datetime.strptime(web[i][0], '%Y-%m-%d %H:%M') for i in range(ct)])
        speed = np.array([float(web[i][1]) for i in range(ct)])
        dir = np.array([float(web[i][2]) for i in range(ct)])

        # unit conversion
        speed = speed/100  # m s-1

        return time, speed, dir

    def get_current(self):
        for i in range(len(self.url)):
            print 'Acquiring data from ' + self.url[i]
            time, speed, dir = self.wget_current(i)
            ctime = nc.date2num(time, 'days since '+self.info['t0'])

            if i==0:
                self.time = time
                self.ctime = ctime
                self.speed = speed
                self.dir = dir
            else:
                if len(ctime)!=len(self.ctime):
                    speed = interp1d(ctime, speed)(self.ctime)
                    dir = interp1d(ctime, dir)(self.ctime)

                self.speed = np.vstack((self.speed, speed))
                self.dir = np.vstack((self.dir, dir))

        return None

    def compute_uv(self):
        self.u = self.speed*np.sin((self.dir-0.)*(np.pi/180))
        self.v = self.speed*np.cos((self.dir-0.)*(np.pi/180))

        self.uraw = self.u.copy()
        self.vraw = self.v.copy()
        self.speedraw = self.speed.copy()
        self.dirraw = self.dir.copy()

        return None

    def filter(self):
        dt = (self.ctime[1]-self.ctime[0])*24.  # hours
        samplefreq = 24./dt  # rad per day
        stopfreq = 24./self.info['Wp_hrs']
        Wp = stopfreq*2./samplefreq
        Ws = 2.*Wp

        n, Wn = buttord(Wp, Ws, 3, 60)
        b, a = butter(n, Wn)
        self.u = filtfilt(b, a, self.u)
        self.v = filtfilt(b, a, self.v)
        U = self.u + 1j*self.v
        self.speed = np.abs(U)
        self.dir = 90-np.angle(U)*180/np.pi

        return None

    def save_data(self):
        fh = nc.Dataset(self.info['filename'], 'w')

        fh.title = 'NOAA ADCP data collection'
        fh.station = self.info['stn']
        fh.station_name = self.info['stn_name']
        fh.lat = self.info['lat']
        fh.lon = self.info['lon']
        fh.deploy_time = self.info['deployed']
        fh.retrieve_time = self.info['retrieved']
        fh.tz = self.info['tz']
        fh.tz_offset = self.info['tz_offset']
        fh.station_url = self.info['stn_url']
        fh.dply_url = self.info['dply_url']
        fh.bins_url = self.info['bins_url']
        fh.bdate = self.info['bdate']
        fh.edate = self.info['edate']
        fh.depth = self.info['depth']
        fh.product = self.info['prod']
        fh.units = self.info['units']
        fh.t0 = self.info['t0']
        fh.bins = self.info['zlevs']

        fh.createDimension('z', self.info['zlevs'])
        fh.createDimension('t')

        fh.createVariable('z', 'd', ('z'))
        fh.createVariable('time', 'd', ('t'))
        fh.createVariable('speedraw', 'd', ('z', 't'))
        fh.createVariable('dirraw', 'd', ('z', 't'))
        fh.createVariable('uraw', 'd', ('z', 't'))
        fh.createVariable('vraw', 'd', ('z', 't'))
        fh.createVariable('speed', 'd', ('z', 't'))
        fh.createVariable('dir', 'd', ('z', 't'))
        fh.createVariable('u', 'd', ('z', 't'))
        fh.createVariable('v', 'd', ('z', 't'))

        fh.variables['z'][:] = self.z
        fh.variables['time'][:] = self.ctime
        fh.variables['speedraw'][:] = self.speedraw
        fh.variables['dirraw'][:] = self.dirraw
        fh.variables['uraw'][:] = self.uraw
        fh.variables['vraw'][:] = self.vraw
        fh.variables['speed'][:] = self.speed
        fh.variables['dir'][:] = self.dir
        fh.variables['u'][:] = self.u
        fh.variables['v'][:] = self.v
        fh.close()

        return None

    def load_data(self):
        fh = nc.Dataset(self.info['filename'], 'r')
        self.z = fh.variables['z'][:]
        self.ctime = fh.variables['time'][:]
        self.speed = fh.variables['speedraw'][:]
        self.dir = fh.variables['dirraw'][:]
        fh.close()

        self.time = nc.num2date(self.ctime, 'days since '+self.info['t0'])

        return None