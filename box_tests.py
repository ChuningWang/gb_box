from nose.tools import *
import numpy as np
import pdb
import netCDF4 as nc

# from box_gb.box_ctd import box_avgbox
# mt, s_box = box_avgbox()

fh = nc.Dataset('./data/s_box.nc')
mt = fh.variables['t'][:]
s_box = fh.variables['s'][:]

testbox = 1
testctd = 0

if testbox == 1:

    from box_gb.box_params import t, F, rtdir, icdir, Cbrs, S0, Sinitbar
    from box_gb import box
    import netCDF4 as nc
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime

    from box_gb.box_params import get_wind
    awnd = get_wind(t)

    S_init = Sinitbar/S0

    S = box.rk_solver(S_init, t, box.box_ode)
    S = S.T*S0
    # pdb.set_trace()

    figs, ax = plt.subplots(4, sharex=True, sharey=True)
    for i in range(4):
        ax[i].plot(mdates.num2date(t), S[:, i*3], 'b')
        ax[i].plot(mdates.num2date(t), S[:, i*3+1], 'g')
        ax[i].plot(mdates.num2date(t), S[:, i*3+2], 'r')
        ax[i].plot(mdates.num2date(mt), s_box[0, i, :], 'ob', ms=3)
        ax[i].plot(mdates.num2date(mt), s_box[1, i, :], 'og', ms=3)
        ax[i].plot(mdates.num2date(mt), s_box[2, i, :], 'or', ms=3)
        ax[i].text(datetime(2006, 01, 01), 26.5, 'Box'+"%1d"%i)

        ax2 = ax[i].twinx()
        # ax2.plot(t, F[:, i], 'grey', lw=.5)
        ax2.plot(mdates.num2date(t), awnd, 'grey', lw=.5, alpha=.1)
        # ax2.set_ylim(0, 3000)
        # ax2.set_yticks([1000, 2000, 3000])

    # ax[3].legend(('Su', 'Si', 'Sd'))
    ax[3].set_xlim(t[0], t[-1])
    ax[3].set_ylim(26, 32)
    ax[3].set_xlabel('Time')
    ax[3].set_ylabel('                                       Salinity [PSU]')
    # ax2.set_ylabel(r'                                        Freshwater Discharge [m$^3\cdot$s$^{-1}$]', color='grey')
    ax[3].set_yticks([27, 29, 31])

    figs.subplots_adjust(hspace=0)

    plt.savefig(rtdir+'S.eps', format='eps')
    plt.close()

if testctd == 1:

    from box_gb.box_params import hu, hi, hd
    from box_gb.box_ctd import box_avgbox

    S_ctd = box_avgbox()

