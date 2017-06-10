from nose.tools import *
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from box_gb.box_const import t, rtdir, icdir, Cbrs
from box_gb import box
import pdb

S0 = np.ones(12)

S = box.rk_solver(S0, t, box.box_ode)
S = S.T
# pdb.set_trace()

figs, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)
ax0.plot(S[:, :3])
ax1.plot(S[:, 3:6])
ax2.plot(S[:, 6:9])
ax3.plot(S[:, 9:12])

ax3.set_ylim(0.7, 1.3)
ax3.legend(('Su','Si','Sd'))


plt.savefig(rtdir+'S.eps',format='eps')
plt.close()
