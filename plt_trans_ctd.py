
# ----------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from gb_toolbox.gb_ctd import rd_ctd, plttrans_ctd, pltmapview_ctd
import gsw

# ---------------------------------------------------------------
# Find indices for each cruise

ctd = rd_ctd('./data/ctd.nc')

sa = gsw.SA_from_SP(ctd['s'], ctd['p'], -136.8, 58.8)
ct = gsw.CT_from_t(sa, ctd['t'], ctd['p'])
rho = gsw.rho(sa, ct, ctd['p'])

dt = np.diff(np.floor(mdates.date2num(ctd['pyt'])))
dt = np.hstack([0, dt])
k1 = np.squeeze(np.where(dt>=7))
k1 = np.hstack([0, k1])
k2 = np.hstack([k1[1:]-1, np.size(dt)+1])
k3 = np.squeeze(np.where((k2-k1)>15))

# Find bottom topography
topo = np.zeros(np.int(np.max(ctd['stn'])))
dd = ctd['d']
dd.mask = ctd['t'].mask
for i in range(np.int(np.max(ctd['stn']))): 
    topo[i] = np.max(dd[:,ctd['stn']==i])
del dd

var = 'rho'
clim=[1000, 1030]

pltfig = 1
if pltfig==1:
    for i in k3[0:1]:
        if var == 'T':
            data1trip = ctd['t'][:, k1[i]:k2[i]]
        elif var == 'S':
            data1trip = ctd['s'][:, k1[i]:k2[i]]
        elif var == 'O2':
            data1trip = ctd['o2'][:, k1[i]:k2[i]]
        elif var == 'Chl':
            data1trip = ctd['f'][:, k1[i]:k2[i]]
        elif var == 'rho':
            data1trip = rho[:, k1[i]:k2[i]]
        stn1trip = ctd['stn'][k1[i]:k2[i]]

        data = np.empty([data1trip.shape[0],ctd['lat_stn'].size])
        data[:] = np.nan
        for j in range(ctd['lat_stn'].size):
            tmp = np.nanmean(data1trip[:,np.where(stn1trip==j)],axis=1) 
            if np.min(tmp.shape)==0:
                data[:,j] = np.nan
            else:
                data[:,j] = np.nanmean(tmp,axis=1)
        del tmp, data1trip, stn1trip
        data = np.ma.masked_invalid(data)

        dis1, dis2, stn_tr1, stn_tr2, ax1, ax2 = plttrans_ctd(
            ctd['lat_stn'],ctd['lon_stn'],np.arange(1, 451),data,clim=clim)
        # Plot topography
        ax1.fill_between(np.concatenate(([0], (dis1[:-1]+dis1[1:])*0.5, [120])),
                         np.concatenate(([topo[stn_tr1[0]]],
                                         topo[stn_tr1[:-1]],
                                         [topo[stn_tr1[-1]]]
                                        )
                                       ),
                         500,
                         facecolor=[.5, .5, .5],edgecolor='none')
        ax2.fill_between(np.concatenate(([0], (dis2[:-1]+dis2[1:])*0.5, [120])),
                         np.concatenate(([topo[stn_tr2[0]]],
                                         topo[stn_tr2[:-1]],
                                         [topo[stn_tr2[-1]]]
                                        )
                                       ),
                         500,
                         facecolor=[.5, .5, .5],edgecolor='none')
        # Set background grey
        ax1.set_axis_bgcolor((.85, .85, .85))
        ax2.set_axis_bgcolor((.85, .85, .85))

        ax1.set_title(var+'__'
                      +"%04d"%ctd['mt_vec'][0,k1[i]]+'_'
                      +"%02d"%ctd['mt_vec'][1,k1[i]]+'_'
                      +"%02d"%ctd['mt_vec'][2,k1[i]],
                      y=1.1)
        plt.savefig('./figs/'+
                    var+'_'+
                    "%04d"%ctd['mt_vec'][0,k1[i]]+
                    "%02d"%ctd['mt_vec'][1,k1[i]]+
                    "%02d"%ctd['mt_vec'][2,k1[i]]+'.eps',format='eps')
        plt.close()
        
pltfig2 = 0
var='Chl'
clevs = np.linspace(0,40,11)
if pltfig2==1:
    for i in k3:
        if var == 'T':
            data1trip = ctd['t'][:,k1[i]:k2[i]]
        elif var == 'S':
            data1trip = ctd['s'][:,k1[i]:k2[i]]
        elif var == 'O2':
            data1trip = ctd['o2'][:,k1[i]:k2[i]]
        elif var == 'Chl':
            data1trip = ctd['f'][:,k1[i]:k2[i]]
        elif var == 'rho':
            data1trip = rho[:,k1[i]:k2[i]]
        stn1trip = ctd['stn'][k1[i]:k2[i]]

        data = np.empty([data1trip.shape[0],ctd['lat_stn'].size])
        data[:] = np.nan
        for j in range(ctd['lat_stn'].size):
            tmp = np.nanmean(data1trip[:,np.where(stn1trip==j)],axis=1) 
            if np.min(tmp.shape)==0:
                data[:,j] = np.nan
            else:
                data[:,j] = np.nanmean(tmp,axis=1)
        del tmp, data1trip, stn1trip
        data = np.ma.masked_invalid(data)
        data[data>np.max(clevs)] = np.max(clevs)
        data[data<np.min(clevs)] = np.min(clevs)
        # Plot map view contour
        k = pltmapview_ctd(ctd['lat_stn'], ctd['lon_stn'], data,
                           d_avg=5, clevs=clevs, cborient='vertical')
        if k:
            plt.title(var+'_'+
                      "%04d"%ctd['mt_vec'][0,k1[i]]+
                      "%02d"%ctd['mt_vec'][1,k1[i]]+
                      "%02d"%ctd['mt_vec'][2,k1[i]]
                     )
            plt.savefig('./figs/mapview/'+
                        var+'_'+
                        "%04d"%ctd['mt_vec'][0,k1[i]]+
                        "%02d"%ctd['mt_vec'][1,k1[i]]+
                        "%02d"%ctd['mt_vec'][2,k1[i]]+'.eps',format='eps')
            plt.close()
 
