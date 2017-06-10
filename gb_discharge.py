import numpy as np

def get_avgbox(boxMethod=1):
    if boxNum==1:
        box0 = np.array([[  -137.7,     59.1  ],
                         [  -137.7,     59.25 ],
                         [  -136.6,     59.25 ],
                         [  -136.15,    58.725],
                         [  -136.65,    58.725],
                         [  -137.15,    58.65 ]])

        box1 = np.array([[  -136.6,     59.25 ],
                         [  -136.3,     59.25 ],
                         [  -135.7,     59.0  ],
                         [  -135.7,     58.725],
                         [  -136.,      58.825],
                         [  -136.15,    58.725],
                         [  -136.3,     58.9  ]])

        box2 = np.array([[  -136.15,    58.725],
                         [  -136.,      58.825],
                         [  -135.7,     58.725],
                         [  -135.6,     58.425],
                         [  -135.975,   58.375],
                         [  -136.0,     58.575]])

        box3 = np.array([[  -136.65,    58.725],
                         [  -136.15,    58.725],
                         [  -136.0,     58.575],
                         [  -135.975,   58.375],
                         [  -136.65,    58.575]])

        box = {'box0': box0,
               'box1': box1,
               'box2': box2,
               'box3': box3
              }

    elif boxNum==2:
        box0 = np.array([[  -137.7,     59.1  ],
                         [  -137.7,     59.25 ],
                         [  -136.6,     59.25 ],
                         [  -136.15,    58.725],
                         [  -136.65,    58.725],
                         [  -137.15,    58.65 ]])

        box1 = np.array([[  -136.6,     59.25 ],
                         [  -136.3,     59.25 ],
                         [  -135.7,     59.0  ],
                         [  -135.7,     58.725],
                         [  -136.,      58.825],
                         [  -136.15,    58.725],
                         [  -136.3,     58.9  ]])

        box2 = np.array([[  -136.15,    58.725],
                         [  -136.,      58.825],
                         [  -135.7,     58.725],
                         [  -135.65,    58.55 ],
                         [  -136.65,    58.55 ],
                         [  -136.65,    58.725]])

        box3 = np.array([[  -136.65,    58.55 ],
                         [  -135.65,    58.55 ],
                         [  -135.6,     58.425],
                         [  -136.0,     58.375]])

        box = {'box0': box0,
               'box1': box1,
               'box2': box2,
               'box3': box3
              }

    return box

#     # ------------------------------------------------------------------------------------------------------
#     # Load discharge data
#     fh = nc.Dataset(pth,'r')
# 
#     t = fh.variables['t'][:]
#     lat = fh.variables['lat'][:]
#     lon = fh.variables['lon'][:]
# 
#     # Get points in boxes
#     hydro_box = np.zeros(lon.shape)
#     p1 = path.Path(box1)
#     p2 = path.Path(box2)
#     p3 = path.Path(box3)
#     p4 = path.Path(box4)
# 
#     for i in range(lon.shape[0]):
#         for j in range(lon.shape[1]):
#             if p1.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 1
#             elif p2.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 2
#             elif p3.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 3
#             elif p4.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 4


