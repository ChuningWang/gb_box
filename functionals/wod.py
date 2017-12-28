""" read in WOD CTD data file. parse data in the Icy Strait. """

import numpy as np
from wodpy import wod
from scipy.interpolate import interp1d
import netCDF4 as nc
from matplotlib import path

# geo range of Icy Strait
c0 = [(-136.10, 58.00),
      (-136.10, 58.45),
      (-136.00, 58.39),
      (-137.00, 58.30),
      (-137.00, 58.00)]

# c0 = [(-136.65, 58.55),
#       (-136.60, 58.35),
#       (-136.50, 57.95),
#       (-135.20, 57.95),
#       (-135.30, 58.50),
#       (-135.60, 58.55),
#       (-135.60, 58.50),
#       (-136.00, 58.35)]

p = path.Path(c0)

depth_sep = 15

time = []
yearday = []
lat = []
lon = []
salt = []
temp = []
zlev = 500
z = np.arange(zlev)

fid = open('/Users/CnWang/Documents/gb_roms/CTDS7513')

n = 0  # counter
t = True
while t:
    pr = wod.WodProfile(fid)
    lon0, lat0 = pr.longitude(), pr.latitude()
    if p.contains_point((lon0, lat0)):
        zpr = pr.z()
        if len(zpr)>1:
            zmin, zmax = zpr[0], zpr[-1]
            zmsk = (z>=zmin) & (z<=zmax)
            spr = np.zeros(zlev)*np.NaN
            spr[zmsk] = interp1d(pr.z(), pr.s())(z[zmsk])
            tpr = np.zeros(zlev)*np.NaN
            tpr[zmsk] = interp1d(pr.z(), pr.t())(z[zmsk])
            if np.nanmax(spr) > 5:
                salt.append(spr)
                temp.append(tpr)
                time.append(nc.date2num(pr.datetime(), 'days since 1900-01-01'))
                yearday.append(pr.datetime().timetuple().tm_yday)
                lat.append(lat0)
                lon.append(lon0)
    t = not pr.is_last_profile_in_file(fid)
    n += 1

salt = np.array(salt)
temp = np.array(temp)
lon = np.array(lon)
lat = np.array(lat)
time = np.array(time)
yearday = np.array(yearday)

# compute seasonal climatology
t = np.array([46, 138, 229, 321])
season_msk = np.zeros((len(time), 4))
season_msk[:, 0] = yearday<=92
season_msk[:, 1] = (yearday>92) & (yearday<=183)
season_msk[:, 2] = (yearday>183) & (yearday<=275)
season_msk[:, 3] = yearday>275

salt_season = np.zeros((4, zlev))
temp_season = np.zeros((4, zlev))

for i in range(4):
    salt_season[i, :] = np.nanmean(salt[season_msk[:, i]==1, :], axis=0)
    temp_season[i, :] = np.nanmean(temp[season_msk[:, i]==1, :], axis=0)

# compute box mean salinity/temperature
Su2 = np.nanmean(salt_season[:, :depth_sep], axis=1)
Si2 = np.nanmean(salt_season[:, depth_sep:], axis=1)
Tu2 = np.nanmean(temp_season[:, :depth_sep], axis=1)
Ti2 = np.nanmean(temp_season[:, depth_sep:], axis=1)

# repeat it three times (last, current, next year)
t = np.concatenate((t-365, t, t+365))
Su2 = np.concatenate((Su2, Su2, Su2))
Si2 = np.concatenate((Si2, Si2, Si2))
Tu2 = np.concatenate((Tu2, Tu2, Tu2))
Ti2 = np.concatenate((Ti2, Ti2, Ti2))

# save the data to txt file
fh = open('../data/gb_icy_strait.txt', 'w')
fh.write('time, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n' % tuple(t))
fh.write('Su2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n' % tuple(Su2))
fh.write('Si2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n' % tuple(Si2))
fh.write('Tu2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n' % tuple(Tu2))
fh.write('Ti2, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n' % tuple(Ti2))
# fh.write('time, Su2, Si2, Tu2, Ti2 \n')
# for i, time in enumerate(t):
#     fh.write('%f, %f, %f, %f, %f \n' % (float(time), Su2[i], Si2[i], Tu2[i], Ti2[i]))
fh.close()
