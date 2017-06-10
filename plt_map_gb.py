# Plot Glacier Bay map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# ------------------------------------------------------------------------------------------------------
# clevs_land = 4*10**np.linspace(0, 3, 101)
# clevs_land = np.linspace(0, 2000, 11)
# clevs_land = []
# clevs_sea = np.linspace(-200, 0, 11)
# clevs_sea = []
# clevs = np.concatenate((clevs_sea, clevs_land))
clevs = np.linspace(-200, 200, 21)
# ------------------------------------------------------------------------------------------------------
ctd = rd_ctd('data/ctd.nc')

# Read topography
fh = nc.Dataset('data/ARDEMv2.0.nc',mode='r')
lon = fh.variables['lon'][:]
lon = lon-360
lat = fh.variables['lat'][:]
z = fh.variables['z'][:]
fh.close()

stn = np.arange(24)

plt.close()
fig = plt.figure()
m = Basemap(llcrnrlon=-138.,llcrnrlat=58,urcrnrlon=-134.5,urcrnrlat=59.4,
            projection='stere',lat_0=58,lon_0=-137.5,
            resolution='i'
           )

xstn, ystn = m(ctd['lon_stn'],ctd['lat_stn'])
m.plot(xstn,ystn,'ok',ms=5)
for i in stn:
    plt.text(xstn[i]+1000,ystn[i],"%02d"%i,fontsize=12,va='center')
# m.drawcoastlines(linewidth=.2)
# m.drawrivers()
# m.fillcontinents(color='lightgrey')
mr = m.drawmeridians(np.arange(-138.,-134.5,0.05),labels=[0,0,0,1],fontsize=6, linewidth=.2)
pr = m.drawparallels(np.arange(58,59.4,0.025),labels=[1,0,0,0],fontsize=6, linewidth=.2)
setlabelrot(mr,-90)

# plot topography

# Use ARDEMv2.0.nc
msk1 = (lon>-138) & (lon<-135)
msk2 = (lat>57.5) & (lat<59.5)

lon2 = lon[msk1]
lat2 = lat[msk2]
lon2, lat2 = np.meshgrid(lon2, lat2)
d = z[msk2,:][:,msk1]
x, y = m(lon2, lat2)
# plot coastline
m.contour(x, y, d, [0], colors='k', linewidths=.6)
dd = m.contour(x, y, d, [-50], colors='k', linewidths=.3, linestyles='solid')
# plt.clabel(dd, fontsize=5)

d[d>0] = d[d>0]/10
d[d<clevs[0]] = clevs[0]
d[d>clevs[-1]] = clevs[-1]

m.contourf(x, y, d, clevs, cmap=plt.cm.coolwarm)
plt.clim(clevs[0], clevs[-1])
cb = plt.colorbar()
cb.set_ticks(clevs)
cb.set_ticklabels(['200', '180', '160', '140', '120', '100', '80', '60', '40', '20', '0', '200', '400', '600', '800', '1000', '1200', '1400', '1600', '1800', '2000'])

# Use GSHHS (Basemap default)

# Plot domain box


# # Plot boxes
# from gb_toolbox.gb_discharge import get_avgbox
# boxes = get_avgbox(boxMethod=1)
# box0 = boxes['box0']
# box1 = boxes['box1']
# box2 = boxes['box2']
# box3 = boxes['box3']
# 
# x, y = m(box0[:,0], box0[:,1])
# m.plot(x, y , 'k', lw=1)
# m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
# x, y = m(box1[:,0], box1[:,1])
# m.plot(x, y , 'k', lw=1)
# m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
# x, y = m(box2[:,0], box2[:,1])
# m.plot(x, y , 'k', lw=1)
# m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
# x, y = m(box3[:,0], box3[:,1])
# m.plot(x, y , 'k', lw=1)
# m.plot([x[-1], x[0]], [y[-1], y[0]], 'k', lw=1)
# 
# x1, y1 = m(-137.6, 59.2)
# plt.text(x1, y1, 'BOX0', fontsize=10)
# x1, y1 = m(-136.5, 59.2)
# plt.text(x1, y1, 'BOX1', fontsize=10)
# x1, y1 = m(-136.6, 58.675)
# plt.text(x1, y1, 'BOX2', fontsize=10)
# x1, y1 = m(-136.4, 58.5)
# plt.text(x1, y1, 'BOX3', fontsize=10)
# 
# 
# # text
# # m_text(-136.2,58.73,'Glacier Bay','rotation',-50,'fontweight','bold')
# # m_text(-136.73,58.1,'Cross Sound','rotation',40,'fontweight','bold')
# # m_text(-136.9,58.75,{'Glacier Bay'; 'National Park'; 'and Preserve'},...
# #               'rotation',-60,'fontweight','bold')
# # m_text(-137.3,59.17,{'A  L  A  S  K  A'},'fontsize',20,'fontweight','bold')
# # m_text(lon_stn+0.03,lat_stn,num2str(transpose(1:24)),'color','r','fontsize',10,'fontweight','bold')
# 
# # x, y = m(-136.2,58.73)
# # plt.text(x,y,'Glacier Bay',rotation=-70,fontweight='bold')

plt.savefig('figs/map_gb_nobox.eps',format='eps')
plt.close()


# # Grid data
# lon_stn = ctd['lon_stn'][stn]
# lat_stn = ctd['lat_stn'][stn]
# 
# numcols, numrows = 300, 300
# xi = np.linspace(lon_stn.min(), lon_stn.max(), numcols)
# yi = np.linspace(lat_stn.min(), lat_stn.max(), numrows)
# xi, yi = np.meshgrid(xi, yi)
# 
# #-- Interpolate at the points in xi, yi
# # "griddata" expects "raw" numpy arrays, so we'll pass in
# # data.x.values instead of just the pandas series data.x
# zi = griddata(lon_stn, lat_stn, np.squeeze(eof1), xi, yi)
# xi, yi = m(xi,yi)
# 
# clevs = np.linspace(0,0.3,31) 
# ci = m.contourf(xi,yi,zi,clevs)
# fig.colorbar(ci)

# plt.show(block=False)
