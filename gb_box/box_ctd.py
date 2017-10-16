from ocean_toolbox import ctd
import numpy as np
# import pdb

def box_avgbox_clim(var='salt', ctddir='./data/ctd.nc', pltfig=-1, svpth=-1, uselowess=-1, filterwindow=30, useallstns=-1):
    # Calculate box climatology (var and freshwater discharge, var can be any variable from CTD measuremsnts.)

    from box_gb.box_params import hu, hi, hd, deltat, boxMethod, t, F
    from scipy.interpolate import interp1d
    import statsmodels.api as sm  # LOWESS filter
    import matplotlib.dates as mdates

    # uselowess = -1

    # Get CTD data
    ctd = rd_ctd(ctddir)
    mt = ctd['mt']  # MATLAB time
    pyt = ctd['pyt']  # Python datetime
    yday = np.array([pyt[i].timetuple().tm_yday for i in range(pyt.size)])  # Year Day
    stn = ctd['stn']  # Station Number
    s = ctd[var]  # var CTD measuremsnts, here use 's' as variable name since it is first developed for salinity

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
        if useallstns == 1:
            # Use every station in each station
            stn2box[(stn==6) | (stn==7) | (stn==8) | (stn==9) | (stn==10) | (stn==11) | (stn==12) | (stn==21)] = 0
            stn2box[(stn==16) | (stn==17) | (stn==18) | (stn==19) | (stn==20)] = 1
            stn2box[(stn==4) | (stn==5) | (stn==13) | (stn==14) | (stn==15) | (stn==11) | (stn==22) | (stn==23)] = 2
            stn2box[(stn==0) | (stn==1) | (stn==2) | (stn==3) | (stn==24)] = 3

        else:
            # Only use station with monthly measurements
            stn2box[(stn==7) | (stn==12)] = 0
            stn2box[(stn==16)] = 1
            stn2box[(stn==4) | (stn==13)] = 2
            stn2box[(stn==0) | (stn==1) | (stn==2) | (stn==3)] = 3

    else:
        print 'The inverse method only works for boxMehtod==1'

    # Average var for upper, intermediate, and deep layer
    ss = np.zeros((3, mt.size))
    for i in range(3):
        ss[0, :] = np.nanmean(s[:hu, :], axis=0)
        ss[1, :] = np.nanmean(s[hu:hu+hi, :], axis=0)
        ss[2, :] = np.nanmean(s[hu+hi:, :], axis=0)

    if uselowess == 1:

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
                lowess = sm.nonparametric.lowess(ss_i[j, :], yday_i, frac=0.1)
                s_box[j, i, :] = interp1d(lowess[:, 0], lowess[:, 1])(yd)

            # River discharge
            lowess = sm.nonparametric.lowess(ff[:, i], yday_ff, frac=0.03)
            f_box[i, :] = interp1d(lowess[:, 0], lowess[:, 1])(yd)

    else:
        # Use bin-average and filtfilt (recommended)

        from scipy.signal import butter, filtfilt

        yday_ff = np.concatenate((yday_f-366, yday_f, yday_f+366))
        ff = np.concatenate((F, F, F), axis = 0)

        for i in range(4):
            # all var measurements in each box
            yday_i = yday[stn2box==i]
            ss_i = ss[:, stn2box==i]

            yday_i = np.concatenate((yday_i-366, yday_i, yday_i+366))
            ss_i = np.concatenate((ss_i, ss_i, ss_i), axis = -1)

            for j in range(yd.size):
                msk = (yday_i>=yd[j]-15) & (yday_i<yd[j]+15)
                s_box[:, i, j] = np.nanmean(ss_i[:, msk], axis=1) 
                msk2 = (yday_ff>=yd[j]-15) & (yday_ff<yd[j]+15)
                f_box[i, j] = np.nanmean(ff[msk2, i]) 

        fs = 1./deltat  # Sampling frequency [day^-1]
        cutoff = 1./filterwindow  # Cutoff frequency [day^-1]
        b, a = butter(5, cutoff/(0.5*fs))

        yd_i = np.concatenate((yd-366, yd, yd+366))
        s_box_i = np.concatenate((s_box, s_box, s_box), axis=-1)
        f_box_i = np.concatenate((f_box, f_box, f_box), axis=-1)

        s_box_i = filtfilt(b, a, s_box_i, axis=-1)
        f_box_i = filtfilt(b, a, f_box_i, axis=-1)

        s_box[:, :, :] = interp1d(yd_i, s_box_i)(yd)
        f_box[:, :] = interp1d(yd_i, f_box_i)(yd)

    # Make a plot of var in each box
    if pltfig == 1:
        import matplotlib.pyplot as plt
        plt.subplot(411)
        plt.plot(yd, s_box[0, 0, :], '--', lw=3)
        plt.plot(yd, s_box[1, 0, :], '--', lw=3)
        plt.plot(yd, s_box[2, 0, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.xlim(0, 366)
        plt.subplot(412)
        plt.plot(yd, s_box[0, 1, :], '--', lw=3)
        plt.plot(yd, s_box[1, 1, :], '--', lw=3)
        plt.plot(yd, s_box[2, 1, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.xlim(0, 366)
        plt.subplot(413)
        plt.plot(yd, s_box[0, 2, :], '--', lw=3)
        plt.plot(yd, s_box[1, 2, :], '--', lw=3)
        plt.plot(yd, s_box[2, 2, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.xlim(0, 366)
        plt.subplot(414)
        plt.plot(yd, s_box[0, 3, :], '--', lw=3)
        plt.plot(yd, s_box[1, 3, :], '--', lw=3)
        plt.plot(yd, s_box[2, 3, :], '--', lw=3)
        plt.ylim(26, 32)
        plt.xlim(0, 366)
        plt.show(block=False)
        plt.savefig(var+'_box.eps', format='eps')

    # Save data
    if svpth != -1:
        import netCDF4 as nc
        # Write data into netCDF file
        print 'Saving data as netCDF4 file in '+svpth+'...'
        f = nc.Dataset(svpth+var+'_box.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay '+var+' & freshwater discharge climatology avaerged with respect to boxes (boxMethod=1)'

        f.createDimension('layer', 3)
        f.createDimension('box', 4)
        f.createDimension('time', None)
        
        t_nc = f.createVariable('t', 'f8', ('time'))
        s_nc = f.createVariable(var+'_box', 'f8', ('layer', 'box', 'time'))
        f_nc = f.createVariable('f_box', 'f8', ('box', 'time'))

        t_nc[:] = yd
        s_nc[:, :, :] = s_box
        f_nc[:, :] = f_box

        f.close()

    return yd, s_box, f_box

def box_avgbox(ctddir='./data/ctd.nc', svpth=-1):

    from box_gb.box_params import hu, hi, hd, icdir, boxMethod
    from scipy import polyfit, polyval

    # Load CTD data
    ctd = rd_ctd(ctddir)
    mt = ctd['mt']
    stn = ctd['stn']
    s = ctd['s']

    st, ed, ltr = get_cruise(ctd)
    mt = 0.5*(mt[st]+mt[ed-2])

    # msk = mt<mdates.date2num(datetime(2009, 01, 01))
    # st = st[msk]
    # ed = ed[msk]
    # mt = mt[msk]

    # Divide into 3 layers, lists them by stations
    ss = np.zeros((3, 25, mt.size))
    ss[:, :, :] = np.NaN

    for i in range(mt.size):
        # Find data from 1 cruise
        s2 = s[:, st[i]:ed[i]]
        stn2 = stn[st[i]:ed[i]]
        # pdb.set_trace()

        for j in range(25):
            msk = stn2==j
            ss[0, j, i] = np.nanmean(s2[:hu, msk])
            ss[1, j, i] = np.nanmean(s2[hu:hu+hi, msk])
            ss[2, j, i] = np.nanmean(s2[hu+hi:, msk])

    s_box = np.zeros((3, 4, mt.size))
    s_box[:, 0, :] = np.squeeze(np.mean(ss[:, [6, 7, 8, 9], :], axis=1))
    s_box[:, 1, :] = np.squeeze(np.mean(ss[:, [16, 17, 18], :], axis=1))
    s_box[:, 2, :] = np.squeeze(np.mean(ss[:, [4, 5, 13], :], axis=1))
    s_box[:, 3, :] = np.squeeze(np.mean(ss[:, [1, 2, 3], :], axis=1))

    # import matplotlib.pyplot as plt
    # plt.plot(mt, ss[1, 1, :], '-o')
    # plt.plot(mt, s_box[1, 3, :], 'o')
    # plt.plot(ss[0, 1, :], s_box[0 ,3, :], 'o')
    # plt.plot(ss[1, 1, :], s_box[1 ,3, :], 'o')
    # plt.show()

    ar = np.zeros((3, 25))
    br = np.zeros((3, 25))
    # stn2box = np.array([3, 3, 3, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 1, 1, 1, 1, 0, 2, 2, 3])
    stn2box = np.array([-1, 3, -1, -1, 2, -1, -1, 0, -1, -1, -1, -1, -1, 2, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1])
    ss2 = np.zeros((3, 25, mt.size))

    for i in range(3):
        for j in range(25):
            idx = ~np.isnan(ss[i, j, :]) & ~np.isnan(s_box[i, stn2box[j], :])
            ar[i, j], br[i, j] = polyfit(np.squeeze(ss[i, j, idx]), np.squeeze(s_box[i, stn2box[j], idx]), 1)
            ss2[i, j, :] = polyval([ar[i, j], br[i, j]], ss[i, j, :])

    for i in range(3):
        for j in range(4):
            for k in range(mt.size):
                kk = np.sum(~np.isnan(ss2[i, stn2box==j, k]))
                if kk>0:
                    s_box[i, j, k] = np.squeeze(np.nanmean(ss2[i, stn2box==j, k]))
                else:
                    s_box[i, j, k] = np.NaN
        
    # For box 3, use intermediate data to substitute deep data, since deep data is too sparse
    s_box[2, 3, :] = s_box[1, 3, :]

    # Remove NaNs
    msk = np.any(np.any(np.isnan(s_box), axis=1), axis=0)
    # pdb.set_trace()
    mt = mt[~msk]
    s_box = s_box[:, :, ~msk]

    # svdata = 1
    if svpth != -1:
        import netCDF4 as nc
        # write data into netCDF file
        print 'Saving data as netCDF4 file in '+svpth+'...'
        f = nc.Dataset(svpth+'s_box.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay salinity avaerged with respect to boxes (boxMethod=2)'

        f.createDimension('layer', 3)
        f.createDimension('box', 4)
        f.createDimension('time', None)
        
        t_nc = f.createVariable('t', 'f8', ('time'))
        s_nc = f.createVariable('s', 'f8', ('layer', 'box', 'time'))

        t_nc[:] = mt
        s_nc[:, :, :] = s_box

        f.close()

    return mt, s_box


