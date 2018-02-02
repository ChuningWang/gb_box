""" Plot Glacier Bay map """

import sys

import numpy as np
from matplotlib import path

from mayavi import mlab
import pyroms
from cmocean import cm

# ------------------- functionals ----------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# ------------------- data prep ------------------------------------------
import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv) > 1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

boxes = {}
boxes['box0'] = np.array([
    [-136.65, 58.65],
    [-137.30, 58.65],
    [-137.30, 59.15],
    [-136.60, 59.15],
    [-136.30, 58.95],
    [-136.20, 58.75]])

boxes['box1'] = np.array([
    [-136.20, 58.75],
    [-136.30, 58.95],
    [-136.60, 59.15],
    [-135.70, 59.15],
    [-135.60, 58.75]])

boxes['box2'] = np.array([
    [-136.65, 58.55],
    [-136.65, 58.65],
    [-136.20, 58.75],
    [-135.60, 58.75],
    [-135.60, 58.55]])

boxes['box3'] = np.array([
    [-136.65, 58.55],
    [-135.60, 58.55],
    [-135.60, 58.50],
    [-136.00, 58.35]])

boxes['box4'] = np.array([
    [-136.65, 58.55],
    [-136.75, 58.25],
    [-136.50, 57.95],
    [-135.20, 57.95],
    [-135.40, 58.55],
    [-135.60, 58.55],
    [-135.60, 58.50],
    [-136.00, 58.35]])

# Read grid
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
z = grd.vgrid.hraw.squeeze()
msk = grd.hgrid.mask_rho
eta, xi = msk.shape

z[z < 0] = 0
z = -z
# z = z/10000
zlevs = 500
z2 = np.arange(0, -1*zlevs, -1)

coords = []
for i in range(eta):
    for j in range(xi):
        coords.append((lon[i, j], lat[i, j]))

b2d = -1 * np.ones((eta, xi))

box_list = boxes.keys()
for box in box_list:
    idx = int(box[-1])
    pp = path.Path(boxes[box])
    box_i = np.reshape(pp.contains_points(coords), (eta, xi))
    b2d[box_i] = idx

b2d[msk == 0] = -1
b3d = np.tile(np.expand_dims(b2d, axis=-1), (1, 1, zlevs))

for i in range(eta):
    for j in range(xi):
        msk = z2 < z[i, j]
        b3d[i, j, msk] = -1

lat3d = np.tile(np.expand_dims(lat, axis=-1), (1, 1, zlevs))
lon3d = np.tile(np.expand_dims(lon, axis=-1), (1, 1, zlevs))
z3d = np.tile(z2, (eta, xi, 1))

# ------------------- make plots -----------------------------------------
lat_min = 58.00
lat_max = 59.20

lon_min = -137.5
lon_max = -135.0

zmin = -450/10000
zmax = 0

mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
mlab.view(azimuth=-90, elevation=18, distance=1)

obj = mlab.contour3d(lon3d, lat3d, z3d/5e3, b3d+1, vmin=-1, vmax=6,
                     colormap='jet', opacity=0.2, transparent=True)
mlab.outline(obj, color=(.7, .7, .7), line_width=1.5)


# # sur = mlab.surf(lon, lat, z,
# #                 extent=[lon_min, lon_max, lat_min, lat_max, zmin, zmax],
# #                 colormap='gray', vmin=zmin, vmax=zmax)
# 
# # mlab.outline(sur, color=(.7, .7, .7), line_width=1.5)
# 
# # mlab.savefig('../figs/map3d.png')
# # mlab.close()
