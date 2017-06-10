# Read freshwater runoff data (hydrology)
# data are stored in /Volumes/R1/ROMS/hydrology/GOA/

import netCDF4 as nc
import glob
from datetime import datetime, timedelta
import pdb
import matplotlib.dates as mdates
import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from gb_ctd import rd_ctd, get_cruise
from matplotlib import path

def get_discharge(pth, savepth, latlim, lonlim):
    # Reading lat lon
    print 'Loading Lat & lon...'
    fh = nc.Dataset(pth+'lat_lon.nc')
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    lon = lon-360

    # Cut out useful portion
    msk_gb = (lat>latlim[0]) & (lat<latlim[1]) & (lon>lonlim[0]) & (lon<lonlim[1])
    msk_lat = np.any(msk_gb,axis=1)
    msk_lon = np.any(msk_gb,axis=0)

    lat = lat[msk_lat, :][:, msk_lon]
    lon = lon[msk_lat, :][:, msk_lon]

    # Load data
    print 'Getting file names...'
    flist = glob.glob(pth+'discharge'+'*.nc')

    rddata = 1
    if rddata == 1:
        # Initiate
        print 'Loading '+flist[0]+'...'
        fh = nc.Dataset(flist[0],mode='r')
        d = fh.variables['discharge'][:]
        d = d[:, :, msk_lat, :][:, :, :, msk_lon]
        t = fh.variables['time'][:]
        t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
        pyt = np.array([t_ini+timedelta(hours=t[i]) for i in range(t.size)])

        for ff in flist[1:]:
            print 'Loading '+ff+'...'
            fh = nc.Dataset(ff,mode='r')
            d_in = fh.variables['discharge'][:]
            d_in = d_in[:, :, msk_lat, :][:, :, :, msk_lon]
            t = fh.variables['time'][:]
            t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
            pyt_in = np.array([t_ini+timedelta(hours=t[i]) for i in range(t.size)])
            
            d = np.concatenate([d, d_in], axis=0)
            pyt = np.concatenate([pyt, pyt_in])

    # mask out invalid data
    print 'Setting invalid data to NaN...'
    d = np.squeeze(d)
    d[d<-1000] = np.nan

    lon_dim = d.shape[1]
    lat_dim = d.shape[2]

    mt = mdates.date2num(pyt)

    # d1 = np.squeeze(d[0,:,:,:])
    # plt.close()
    # plt.contourf(lon,lat,d1)
    # plt.show(block=False)

    svdata = 1
    if svdata == 1:
        # write data into netCDF file
        print 'Saving data as netCDF4 file...'
        f = nc.Dataset(savepth+'discharge_gb.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay river discharge and deglaciation'

        f.createDimension('time', None)
        f.createDimension('lat', lat_dim)
        f.createDimension('lon', lon_dim)

        t_nc = f.createVariable('t', 'f8', ('time'))
        lat_nc = f.createVariable('lat', 'f8', ('lon', 'lat'))
        lon_nc = f.createVariable('lon', 'f8', ('lon', 'lat'))
        d_nc = f.createVariable('discharge', 'f8', ('time', 'lon', 'lat'))

        t_nc[:] = mt
        lat_nc[:, :] = lat
        lon_nc[:, :] = lon
        d_nc[:, :, :] = d

        f.close()

def get_discharge_avgbox(pth, savepth=-1):
    from gb_discharge import get_avgbox
    boxes = get_avgbox(boxMethod=2)
    box0 = boxes['box0']
    box1 = boxes['box1']
    box2 = boxes['box2']
    box3 = boxes['box3']
    # ------------------------------------------------------------------------------------------------------
    # Load discharge data
    fh = nc.Dataset(pth,'r')

    t = fh.variables['t'][:]
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]

    # Get points in boxes
    hydro_box = np.ones(lon.shape)*(-1)
    p0 = path.Path(box0)
    p1 = path.Path(box1)
    p2 = path.Path(box2)
    p3 = path.Path(box3)

    for i in range(lon.shape[0]):
        for j in range(lon.shape[1]):
            if p0.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 0
            elif p1.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 1
            elif p2.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 2
            elif p3.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 3

    # ------------------------------------------------------------------------------------------------------
    # Divide GB into several hydro regions
    avg_runoff = 1
    if avg_runoff == 1:
        d = np.empty((t.size, 4))
        d[:] = np.NaN
        for i in range(t.size):
            d0 = np.squeeze(fh.variables['discharge'][i, :, :])
            d[i, 0] = np.nanmean(d0[hydro_box==0])
            d[i, 1] = np.nanmean(d0[hydro_box==1])
            d[i, 2] = np.nanmean(d0[hydro_box==2])
            d[i, 3] = np.nanmean(d0[hydro_box==3])

        # plt.close()
        # plt.figure()
        # plt.plot(pyt, d)
        # plt.legend(('1','2','3','4'))
        # plt.savefig('../figs/runoff.eps', format='eps')

        discharge = {'mt':          t,
                     'box0':        box0,
                     'box1':        box1,
                     'box2':        box2,
                     'box3':        box3,
                     'discharge':   d
                    }

    if savepth!=-1:
        f = nc.Dataset(savepth+'discharge_gb_box.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay freshwater discharge and deglaciation, sum of each box'

        f.createDimension('time', None)
        f.createDimension('box', None)

        t_nc = f.createVariable('t', 'f8', ('time'))
        d_nc = f.createVariable('discharge', 'f8', ('time', 'box'))

        t_nc[:] = t
        d_nc[:, :] = d
        t_nc.units = 'days since 0001-01-01'
        d_nc.units = 'm^3s^-1'

        f.close()

    return discharge

def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

pth = '/Volumes/R1/ROMS/hydrology/GOA/'
savepth = '/Users/chuning/projects/glacierbay/python/data/'
latlim = [58.3, 59.3]
lonlim = [-138.0, -135.5]

# ------------------------------------------------------------------------------------------------------
# get_discharge(pth, savepth, latlim, lonlim)

# ------------------------------------------------------------------------------------------------------
discharge = get_discharge_avgbox('data/discharge_gb.nc', savepth)
t = discharge['mt']
box1 = discharge['box1']
box2 = discharge['box2']
box3 = discharge['box3']
box4 = discharge['box4']

# ------------------------------------------------------------------------------------------------------
ctd = rd_ctd('data/ctd.nc')
st, ed, ltr = get_cruise(ctd)
t_ctd = np.floor(mdates.date2num(ctd['pyt']))
t_ctd = t_ctd[st]
lat_stn = ctd['lat_stn']
lon_stn = ctd['lon_stn']

fh = nc.Dataset('data/discharge_gb.nc')
lon = fh.variables['lon'][:]
lat = fh.variables['lat'][:]

# ------------------------------------------------------------------------------------------------------
pltfig_tst = 0
pltlog = 1
if pltfig_tst == 1:
    d0 = np.squeeze(fh.variables['discharge'][0,:,:])
    plt.close()
    fig = plt.figure()
    m = Basemap(llcrnrlon=-138.,llcrnrlat=58,urcrnrlon=-135.5,urcrnrlat=59.4,
                projection='stere',lat_0=58,lon_0=-137.5,
                resolution='i'
               )

    # pdb.set_trace()
    x, y = m(lon_stn,lat_stn)
    m.plot(x,y,'ok',markersize=3)
    m.drawcoastlines(linewidth=.5)

    mr = m.drawmeridians(np.arange(-138.,-135.5,0.05),labels=[0,0,0,1],fontsize=6, linewidth=.2)
    pr = m.drawparallels(np.arange(58,59.4,0.025),labels=[1,0,0,0],fontsize=6, linewidth=.2)
    setlabelrot(mr,-90)

    # Plot boxes
    x, y = m(box1[:,0], box1[:,1])
    m.plot(x, y , 'k', lw=1)
    m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
    x, y = m(box2[:,0], box2[:,1])
    m.plot(x, y , 'k', lw=1)
    m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
    x, y = m(box3[:,0], box3[:,1])
    m.plot(x, y , 'k', lw=1)
    m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
    x, y = m(box4[:,0], box4[:,1])
    m.plot(x, y , 'k', lw=1)
    m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)

    xi, yi = m(lon, lat)

    if pltlog == 1:
        clevs = np.arange(-1, 8)
        d0 = np.log10(d0)
    else:
        clevs = np.linspace(0, 1, 11)
        d0 = d0/1e7

    d0[d0<clevs[0]] = clevs[0]
    d0[d0>clevs[-1]] = clevs[-1]
    m.contourf(xi, yi, d0, clevs)

    if pltlog == 1:
        cbar = plt.colorbar(label=r'Freshwater Runoff [Log Scale, m$^3\cdot$s$^{-1}$]')
    else:
        cbar = plt.colorbar(label=r'Freshwater Runoff [10$^7$ m$^3\cdot$s$^{-1}$]')

    plt.title('Runoff')
    # plt.show(block=False)
    plt.savefig('../figs/runoff_t0.eps', format='eps')

# ------------------------------------------------------------------------------------------------------
pltfig = 0
pltlog = 1
if pltfig == 1:
    fh = nc.Dataset('data/discharge_gb.nc')
    for tt in t_ctd:
        idx = np.where(t==tt)[0]
        if idx.shape[0]==0:
            print('No runoff data on '+mdates.num2date(tt).strftime('%Y-%m-%d'))
            break
        else:
            print('Found runoff data on '+mdates.num2date(tt).strftime('%Y-%m-%d'))
            d0 = np.squeeze(fh.variables['discharge'][idx, :, :])

        print('Plotting...')
        plt.close()
        fig = plt.figure()
        m = Basemap(llcrnrlon=-138.,llcrnrlat=58,urcrnrlon=-135.5,urcrnrlat=59.4,
                    projection='stere',lat_0=58,lon_0=-137.5,
                    resolution='i'
                   )

        # pdb.set_trace()
        x, y = m(lon_stn,lat_stn)
        m.plot(x,y,'ok',markersize=3)
        m.drawcoastlines(linewidth=.5)

        mr = m.drawmeridians(np.arange(-138.,-135.5,0.05),labels=[0,0,0,1],fontsize=6, linewidth=.2)
        pr = m.drawparallels(np.arange(58,59.4,0.025),labels=[1,0,0,0],fontsize=6, linewidth=.2)
        setlabelrot(mr,-90)

        # Plot boxes
        x, y = m(box1[:,0], box1[:,1])
        m.plot(x, y , 'k', lw=.5)
        m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=.5)
        x, y = m(box2[:,0], box2[:,1])
        m.plot(x, y , 'k', lw=.5)
        m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=.5)
        x, y = m(box3[:,0], box3[:,1])
        m.plot(x, y , 'k', lw=.5)
        m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=.5)
        x, y = m(box4[:,0], box4[:,1])
        m.plot(x, y , 'k', lw=.5)
        m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=.5)

        xi, yi = m(fh.variables['lon'][:], fh.variables['lat'][:])

        if pltlog == 1:
            clevs = np.arange(-1, 8)
            d0 = np.log10(d0)
        else:
            clevs = np.linspace(0, 1, 11)
            d0 = d0/1e7

        d0[d0<clevs[0]] = clevs[0]
        d0[d0>clevs[-1]] = clevs[-1]
        m.contourf(xi, yi, d0, clevs)

        if pltlog == 1:
            cbar = plt.colorbar(label=r'Freshwater Runoff [Log Scale, m$^3\cdot$s$^{-1}$]')
        else:
            cbar = plt.colorbar(label=r'Freshwater Runoff [10$^7$ m$^3\cdot$s$^{-1}$]')

        plt.title('Runoff '+mdates.num2date(tt).strftime('%Y-%m-%d'))
        # plt.show(block=False)
        plt.savefig('../figs/runoff/r'+mdates.num2date(tt).strftime('%Y%m%d')+'.eps', format='eps')


