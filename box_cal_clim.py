# Calculate box climatology (salinity, River discharge, ect)

from box_gb.box_params import hu, hi, hd, deltat, boxMethod, t, F
import numpy as np
from gb_toolbox.gb_ctd import rd_ctd

def box_cal_clim(ctddir='./data/ctd.nc'):
    from scipy.interpolate import interp1d
    import statsmodels.api as sm  # LOWESS filter
    import matplotlib.dates as mdates

    # Get CTD data
    ctd = rd_ctd(ctddir)
    mt = ctd['mt']  # MATLAB time
    pyt = ctd['pyt']  # Python datetime
    yday = np.array([pyt[i].timetuple().tm_yday for i in range(pyt.size)])  # Year Day
    stn = ctd['stn']  # Station Number
    s = ctd['s']  # Salinity

    # Convert Freshwater discharge MATLAB time into Python datetime
    pyt_f = np.array(mdates.num2date(t))
    yday_f = np.array([pyt_f[i].timetuple().tm_yday for i in range(pyt_f.size)])  # Year Day

    # Initiate
    yd = np.arange(0, 366, deltat)  # Time grid
    s_box = np.zeros((3, 4, yd.size))
    f_box = np.zeros((4, yd.size))

    # Find which box each CTD cast belongs to
    # stn2box = np.array([3, 3, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 1, 1, 1, 1, 0, 2, 2, 3])
    stn2box = -1*np.ones(stn.size)  # Initiate
    if boxMethod == 1:
        stn2box[(stn==6) | (stn==7) | (stn==8) | (stn==9) | (stn==10) | (stn==11) | (stn==12) | (stn==21)] = 0
        stn2box[(stn==16) | (stn==17) | (stn==18) | (stn==19) | (stn==20)] = 1
        stn2box[(stn==4) | (stn==5) | (stn==13) | (stn==14) | (stn==15) | (stn==11) | (stn==22) | (stn==23)] = 2
        stn2box[(stn==0) | (stn==1) | (stn==2) | (stn==3) | (stn==24)] = 3
    else:
        print 'The inverse method only works for boxMehtod==1'

    # Average Salinity for upper, intermediate, and deep layer
    ss = np.zeros((3, mt.size))
    for i in range(3):
        ss[0, :] = np.nanmean(s[:hu, :], axis=0)
        ss[1, :] = np.nanmean(s[hu:hu+hi, :], axis=0)
        ss[2, :] = np.nanmean(s[hu+hi:, :], axis=0)


    # Use LOWESS filter to smooth the data and then interpolate onto yd
    yday_ff = np.concatenate((yday_f-366, yday_f, yday_f+366))
    ff = np.concatenate((F, F, F), axis = 0)

    for i in range(4):
        # Salinity
        ss_i = ss[:, stn2box==i]
        yday_i = yday[stn2box==i]

        yday_i = np.concatenate((yday_i-366, yday_i, yday_i+366))
        ss_i = np.concatenate((ss_i, ss_i, ss_i), axis = 1)

        for j in range(3):
            lowess = sm.nonparametric.lowess(ss_i[j, :], yday_i, frac=0.06)
            s_box[j, i, :] = interp1d(lowess[:, 0], lowess[:, 1])(yd)

        # River discharge
        lowess = sm.nonparametric.lowess(ff[:, i], yday_ff, frac=0.03)
        f_box[i, :] = interp1d(lowess[:, 0], lowess[:, 1])(yd)

    # Make a plot of salinities in each box
    pltfig = 0
    if pltfig == 1:
        import matplotlib.pyplot as plt
        plt.subplot(411)
        plt.plot(yd, s_box[0, 0, :], '--', lw=3)
        plt.plot(yd, s_box[1, 0, :], '--', lw=3)
        plt.plot(yd, s_box[2, 0, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.subplot(412)
        plt.plot(yd, s_box[0, 1, :], '--', lw=3)
        plt.plot(yd, s_box[1, 1, :], '--', lw=3)
        plt.plot(yd, s_box[2, 1, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.subplot(413)
        plt.plot(yd, s_box[0, 2, :], '--', lw=3)
        plt.plot(yd, s_box[1, 2, :], '--', lw=3)
        plt.plot(yd, s_box[2, 2, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.subplot(414)
        plt.plot(yd, s_box[0, 3, :], '--', lw=3)
        plt.plot(yd, s_box[1, 3, :], '--', lw=3)
        plt.plot(yd, s_box[2, 3, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.show(block=False)
        plt.savefig('s_box.eps', format='eps')

    # Save data
    svdata = 1
    if svdata == 1:
        import netCDF4 as nc
        # Write data into netCDF file
        print 'Saving data as netCDF4 file...'
        f = nc.Dataset('./data/s_box.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay salinity & freshwater discharge climatology avaerged with respect to boxes (boxMethod=1)'

        f.createDimension('layer', 3)
        f.createDimension('box', 4)
        f.createDimension('time', None)
        
        t_nc = f.createVariable('t', 'f8', ('time'))
        s_nc = f.createVariable('s_box', 'f8', ('layer', 'box', 'time'))
        f_nc = f.createVariable('f_box', 'f8', ('box', 'time'))

        t_nc[:] = yd
        s_nc[:, :, :] = s_box
        f_nc[:, :] = f_box

        f.close()

    return yd, s_box, f_box
