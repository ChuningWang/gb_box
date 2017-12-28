# Inverse box model
# 2016/09/06
# Add a equation for temperature for the bottom box

def box_inverse(cal_clim=-1, svpth=-1, filterwindow=90):

    import numpy as np
    from box_gb.box_params import Aui, Aid, Vu, Vi, Vd, deltaSui, deltaSid, damp_mixing
    from box_gb.box_params import get_pacific
    from scipy.optimize import lsq_linear  # Solve linear system with constrains

    # cal_clim = 0

    # Calculate/Read salinities and freshwater discharge in each box
    # The function box_avgbox_clim works fine on pool-97, but does not work on my laptop.
    # Thus use cal_clim = - 1 unless you want to change the averaging algorithm.
    if cal_clim == 1:
        from box_gb.box_ctd import box_avgbox_clim
        mt, s, F = box_avgbox_clim(svpth='./data/', filterwindow=filterwindow)
    else:
        import netCDF4 as nc
        fh = nc.Dataset('./data/s_box.nc', mode='r')
        mt = fh.variables['t'][:]
        s = fh.variables['s_box'][:]
        F = fh.variables['f_box'][:]

    # ----------------------------------------------------------------------------------------------
    dmt = np.diff(mt)*24*60*60  # [s]
    dsdt = np.diff(s, axis=2)/dmt  # [PSU s^-1]

    mt = 0.5*(mt[:-1]+mt[1:])
    s = 0.5*(s[:, :, :-1]+s[:, :, 1:])  # [PSU]
    F = 0.5*(F[:, :-1]+F[:, 1:])  # [m^3s^-1]

    # Pacific Salinity (sp)
    # sp = get_pacific(mt, lat_sp=56.875, lon_sp=-137.625)
    sp = get_pacific(mt, filterwindow=filterwindow)  # [PSU]

    X = np.zeros((12, mt.size))*np.nan

    for i in range(mt.size):
        # Solve linear system aX=b with constraint X>0
        a = np.zeros((11, 11))

        # --------------------------------------------------------------------------------
        # This version is created on 2016/08/23 -- see my notes
        # 2016/08/24
        # Change this from a 11-box model into a 12-box model -- see my notes
                
        a[0, 0]   = s[1, 0, i]-s[0, 0, i]
        a[0, 4]   = s[1, 0, i]-s[0, 0, i]
        a[1, 0]   = s[2, 0, i]-s[1, 0, i]
        a[1, 4]   = s[0, 0, i]-s[1, 0, i]
        a[1, 5]   = s[2, 0, i]-s[1, 0, i]
        a[2, 0]   = s[2, 2, i]-s[2, 0, i]
        a[2, 5]   = s[1, 0, i]-s[2, 0, i]
        a[3, 1]   = s[1, 1, i]-s[0, 1, i]
        a[3, 6]   = s[1, 1, i]-s[0, 1, i]
        a[4, 1]   = s[2, 1, i]-s[1, 1, i]
        a[4, 6]   = s[0, 1, i]-s[1, 1, i]
        a[4, 7]   = s[2, 1, i]-s[1, 1, i]
        a[5, 1]   = s[2, 2, i]-s[2, 1, i]
        a[5, 7]   = s[1, 1, i]-s[2, 1, i]
        a[6, 0]   = s[0, 0, i]-s[1, 2, i]
        a[6, 1]   = s[0, 1, i]-s[1, 2, i]
        a[6, 2]   = s[1, 2, i]-s[0, 2, i]
        a[6, 8]   = s[1, 2, i]-s[0, 2, i]
        a[7, 0]   = s[1, 2, i]-s[2, 2, i]
        a[7, 1]   = s[1, 2, i]-s[2, 2, i]
        a[7, 2]   = s[2, 2, i]-s[1, 2, i]
        a[7, 8]   = s[0, 2, i]-s[1, 2, i]
        a[7, 9]   = s[2, 2, i]-s[1, 2, i]
        a[8, 2]   = s[1, 3, i]-s[2, 2, i]
        a[8, 9]   = s[1, 2, i]-s[2, 2, i]
        a[9, 2]   = s[0, 2, i]-s[1, 3, i]
        a[9, 3]   = s[1, 3, i]-s[0, 3, i]
        a[9, 10]  = s[1, 3, i]-s[0, 3, i]
        a[10, 3]  = sp[i]-s[1, 3, i]
        a[10, 10] = s[0, 3, i]-s[1, 3, i]

        b = np.array([Vu[0]*dsdt[0, 0, i]+F[0, i]*s[1, 0, i],
                      Vi[0]*dsdt[1, 0, i]+F[0, i]*(s[2, 0, i]-s[1, 0, i]),
                      Vd[0]*dsdt[2, 0, i]+F[0, i]*(s[2, 2, i]-s[2, 0, i]),
                      Vu[1]*dsdt[0, 1, i]+F[1, i]*s[1, 1, i],
                      Vi[1]*dsdt[1, 1, i]+F[1, i]*(s[2, 1, i]-s[1, 1, i]),
                      Vd[1]*dsdt[2, 1, i]+F[1, i]*(s[2, 2, i]-s[2, 0, i]),
                      Vu[2]*dsdt[0, 2, i]+F[2, i]*s[1, 2, i],
                      Vi[2]*dsdt[1, 2, i]+F[2, i]*(s[2, 2, i]-s[1, 2, i]),
                      Vd[2]*dsdt[2, 2, i]+(F[0, i]+F[1, i]+F[2, i])*(s[1, 3, i]-s[2, 2, i]),
                      Vu[3]*dsdt[0, 3, i]+F[3, i]*s[1, 3, i],
                      Vi[3]*dsdt[1, 3, i]+(F[0, i]+F[1, i]+F[2, i]+F[3, i])*(sp[i]-s[1, 3, i]),
                     ])

        # --------------------------------------------------------------------------------
        # This is the original version

        # a[0, 0]   =  (s[0, 2, i]-s[0, 0, i])*(s[1, 0, i]-s[0, 0, i])
        # a[0, 3]   =  Aui[0]*damp_ui[0]*(s[1, 0, i]-s[0, 0, i])
        # a[1, 0]   =  (s[0, 2, i]-s[0, 0, i])*(s[1, 2, i]-s[1, 0, i])
        # a[1, 3]   = -Aui[0]*damp_ui[0]*(s[1, 0, i]-s[0, 0, i])
        # a[1, 4]   =  Aid[0]*damp_id[0]*(s[2, 0, i]-s[1, 0, i])
        # a[2, 4]   = -Aid[0]*damp_id[0]*(s[2, 0, i]-s[1, 0, i])
        # a[3, 1]   =  (s[0, 2, i]-s[0, 1, i])*(s[1, 1, i]-s[1, 0, i])
        # a[3, 5]   =  Aui[1]*damp_ui[1]*(s[1, 1, i]-s[0, 1, i])
        # a[4, 1]   =  (s[0, 2, i]-s[0, 1, i])*(s[1, 2, i]-s[1, 1, i])
        # a[4, 5]   = -Aui[1]*damp_ui[1]*(s[1, 1, i]-s[0, 1, i])
        # a[4, 6]   =  Aid[1]*damp_id[1]*(s[2, 1, i]-s[1, 1, i])
        # a[5, 6]   = -Aid[1]*damp_id[1]*(s[2, 1, i]-s[1, 1, i])
        # a[6, 0]   =  (s[0, 2, i]-s[0, 0, i])*(s[0, 0, i]-s[1, 2, i])
        # a[6, 1]   =  (s[0, 2, i]-s[0, 1, i])*(s[0, 1, i]-s[1, 2, i])
        # a[6, 2]   =  (s[0, 3, i]-s[0, 2, i])*(s[1, 2, i]-s[0, 2, i])
        # a[6, 7]   =  Aui[2]*damp_ui[2]*(s[1, 2, i]-s[0, 2, i])
        # a[7, 2]   =  (s[0, 3, i]-s[0, 2, i])*(s[1, 3, i]-s[1, 2, i])
        # a[7, 7]   = -Aui[2]*damp_ui[2]*(s[1, 2, i]-s[0, 2, i])
        # a[7, 8]   =  Aid[2]*damp_id[2]*(s[2, 2, i]-s[1, 2, i])
        # a[8, 8]   = -Aid[2]*damp_id[2]*(s[2, 2, i]-s[1, 2, i])
        # a[9, 2]   =  (s[0, 3, i]-s[0, 2, i])*(s[0, 2, i]-s[0, 3, i])
        # a[9, 9]   =  Aui[3]*damp_ui[3]*(s[1, 3, i]-s[0, 3, i])
        # a[10, 2]  =  (s[0, 3, i]-s[0, 2, i])*(Spbar-s[1, 3, i])
        # a[10, 9]  = -Aui[3]*damp_ui[3]*(s[1, 3, i]-s[0, 3, i])
        # a[10, 10] =  Aid[3]*damp_id[3]*(s[2, 3, i]-s[1, 3, i])
        # a[10, 11] = Vi[3]*(Spbar-s[1, 3, i])
        # a[11, 10] = -Aid[3]*damp_id[3]*(s[2, 3, i]-s[1, 3, i])

        # b = np.array([Vu[0]*dsdt[0, 0, i]+F[0, i]*s[1, 0, i],
        #               Vi[0]*dsdt[1, 0, i]+F[0, i]*(s[1, 2, i]-s[1, 0 ,i]),
        #               Vd[0]*dsdt[2, 0, i],
        #               Vu[1]*dsdt[0, 1, i]+F[1, i]*s[1, 1, i],
        #               Vi[1]*dsdt[1, 1, i]+F[1, i]*(s[1, 2, i]-s[1, 1, i]),
        #               Vd[1]*dsdt[2, 1, i],
        #               Vu[2]*dsdt[0, 2, i],
        #               Vi[2]*dsdt[1, 2, i]+(F[0, i]+F[1, i]+F[2, i])*(s[1, 3, i]-s[1, 2, i]),
        #               Vd[2]*dsdt[2, 2, i],
        #               Vu[3]*dsdt[0, 3, i]+F[3, i]*s[0, 3, i],
        #               Vi[3]*dsdt[1, 3, i]+(F[0, i]+F[1, i]+F[2, i]+F[3, i])*(Spbar-s[1, 3, i]),
        #               Vd[3]*dsdt[2, 3, i],
        #              ])

        # X[:-1, i] = np.linalg.solve(a, b)

        X[:-1, i] = lsq_linear(a, b, bounds=(np.zeros(11), np.ones(11)*np.inf)).x

    # ----------------------------------------------------------------------------------------------
    # Save data
    # svdata = 1
    if svpth != -1:
        import netCDF4 as nc
        # Write data into netCDF file
        print 'Saving data as netCDF4 file to '+svpth+'...'
        f = nc.Dataset(svpth+'inverse.nc', 'w', format='NETCDF4')
        f.description = 'Outout of the inverse box model (boxMethod=1)'

        f.createDimension('time', None)
        f.createDimension('layer', 3)
        f.createDimension('box', 4)
        
        t_nc = f.createVariable('t', 'f8', ('time'))
        q02_nc = f.createVariable('Q02', 'f8', ('time'))
        q12_nc = f.createVariable('Q12', 'f8', ('time'))
        q23_nc = f.createVariable('Q23', 'f8', ('time'))
        q3p_nc = f.createVariable('Q3p', 'f8', ('time'))
        mix_ui_nc = f.createVariable('Wui', 'f8', ('box', 'time'))
        mix_id_nc = f.createVariable('Wid', 'f8', ('box', 'time'))
        s_nc = f.createVariable('s_box', 'f8', ('layer', 'box', 'time'))
        f_nc = f.createVariable('f_box', 'f8', ('box', 'time'))
        sp_nc = f.createVariable('sp', 'f8', ('time'))

        t_nc[:] = mt
        s_nc[:, :, :] = s
        f_nc[:, :] = F
        sp_nc[:] = sp
        q02_nc[:] = X[0, :]
        q12_nc[:] = X[1, :]
        q23_nc[:] = X[2, :]
        q3p_nc[:] = X[3, :]
        mix_ui_nc[:, :] = X[[4, 6, 8, 10], :]
        mix_id_nc[:, :] = X[[5, 7, 9, 11], :]

        f.close()

    inv = {'mt':    mt,
           's_box': s,
           'f_box': F,
           'sp':    sp,
           'Q02':   X[0, :],
           'Q12':   X[1, :],
           'Q23':   X[2, :],
           'Q3p':   X[3, :],
           'Wui':   X[[4, 6, 8, 10], :],
           'Wid':   X[[5, 7, 9, 11], :]
          }

    return inv


def box_inverse2(cal_clim=-1, svpth=-1, filterwindow=90):

    import numpy as np
    from box_gb.box_params import Aui, Aid, Vu, Vi, Vd, deltaSui, deltaSid, damp_mixing
    from box_gb.box_params import get_pacific
    from scipy.optimize import lsq_linear  # Solve linear system with constrains

    heatflux = 0  # heat budget equation is a failure -- keep heatflux = 0
    tidemixing = 1

    # Calculate/Read salinities and freshwater discharge in each box
    # The function box_avgbox_clim works fine on pool-97, but does not work on my laptop.
    # Thus use cal_clim = - 1 unless you want to change the averaging algorithm.
    if cal_clim == 1:
        from box_gb.box_ctd import box_avgbox_clim
        mt, s, F = box_avgbox_clim(svpth='./data/', filterwindow=filterwindow)
        mt, t, F = box_avgbox_clim(var='t', svpth='./data/', filterwindow=filterwindow)
    else:
        import netCDF4 as nc
        fh = nc.Dataset('./data/s_box.nc', mode='r')
        mt = fh.variables['t'][:]
        s = fh.variables['s_box'][:]
        F = fh.variables['f_box'][:]
        fh.close()
        
        fh = nc.Dataset('./data/t_box.nc', mode='r')
        t = fh.variables['t_box'][:]
        fh.close()

    su = s[0, :, :]
    si0 = (s[1, 0, :]*Vi[0]+s[1, 1, :]*Vi[1]+s[1, 2, :]*Vi[2])/(Vi[0]+Vi[1]+Vi[2])
    si3 = s[1, 3, :]
    sd = (s[2, 0, :]*Vd[0]+s[2, 1, :]*Vd[1]+s[2, 2, :]*Vd[2])/(Vd[0]+Vd[1]+Vd[2])
    Vi0 = Vi[0]+Vi[1]+Vi[2]
    Vi3 = Vi[3]
    Vd0 = Vd[0]+Vd[1]+Vd[2]
    Aid0 = Aid[0]+Aid[1]+Aid[2]
    # ----------------------------------------------------------------------------------------------
    dmt = np.diff(mt)*24*60*60  # [s]
    dsdt = np.diff(s, axis=2)/dmt  # [PSU s^-1]
    dsudt = dsdt[0, :, :]
    dsi0dt = (dsdt[1, 0, :]*Vi[0]+dsdt[1, 1, :]*Vi[1]+dsdt[1, 2, :]*Vi[2])/(Vi[0]+Vi[1]+Vi[2])
    dsi3dt = dsdt[1, 3, :]
    dsddt = (dsdt[2, 0, :]*Vd[0]+dsdt[2, 1, :]*Vd[1]+dsdt[2, 2, :]*Vd[2])/(Vd[0]+Vd[1]+Vd[2])
    # ----------------------------------------------------------------------------------------------
    mt = 0.5*(mt[:-1]+mt[1:])
    su = 0.5*(su[:, :-1]+su[:, 1:])  # [PSU]
    si0 = 0.5*(si0[:-1]+si0[1:])  # [PSU]
    si3 = 0.5*(si3[:-1]+si3[1:])  # [PSU]
    sd = 0.5*(sd[:-1]+sd[1:])  # [PSU]
    s = 0.5*(s[:, :, :-1]+s[:, :, 1:])  # [PSU]
    F = 0.5*(F[:, :-1]+F[:, 1:])  # [m^3s^-1]

    # Pacific Salinity (sp)
    # sp = get_pacific(mt, lat_sp=56.875, lon_sp=-137.625, filterwindow=filterwindow)
    sp = get_pacific(mt, filterwindow=filterwindow)  # [PSU]

    # ----------------------------------------------------------------------------------------------
    if heatflux == 1:
        # Calculate temperature
        tu = t[0, :, :]
        ti0 = (t[1, 0, :]*Vi[0]+t[1, 1, :]*Vi[1]+t[1, 2, :]*Vi[2])/(Vi[0]+Vi[1]+Vi[2])
        ti3 = t[1, 3, :]
        td = (t[2, 0, :]*Vd[0]+t[2, 1, :]*Vd[1]+t[2, 2, :]*Vd[2])/(Vd[0]+Vd[1]+Vd[2])

        dtdt = np.diff(t, axis=2)/dmt  # [PSU s^-1]
        dtudt = dtdt[0, :, :]
        dti0dt = (dtdt[1, 0, :]*Vi[0]+dtdt[1, 1, :]*Vi[1]+dtdt[1, 2, :]*Vi[2])/(Vi[0]+Vi[1]+Vi[2])
        dti3dt = dtdt[1, 3, :]
        dtddt = (dtdt[2, 0, :]*Vd[0]+dtdt[2, 1, :]*Vd[1]+dtdt[2, 2, :]*Vd[2])/(Vd[0]+Vd[1]+Vd[2])

        tu = 0.5*(tu[:, :-1]+tu[:, 1:])  # [PSU]
        ti0 = 0.5*(ti0[:-1]+ti0[1:])  # [PSU]
        ti3 = 0.5*(ti3[:-1]+ti3[1:])  # [PSU]
        td = 0.5*(td[:-1]+td[1:])  # [PSU]
        t = 0.5*(t[:, :, :-1]+t[:, :, 1:])  # [PSU]

        tp = get_pacific(mt, var='t', filterwindow=filterwindow)  # [PSU]

    # ----------------------------------------------------------------------------------------------
    X = np.zeros((7, mt.size))*np.nan

    for i in range(mt.size):
        # Solve linear system aX=b with constraint X>0
        if heatflux == 1:
            a = np.zeros((10, 7))
        else:
            a = np.zeros((7, 7))

        # --------------------------------------------------------------------------------
        # This version is created on 2016/09/01 -- see my notes
        # 2016/09/01
        # Combine intermediate and deep layer for box 0, 1, and 2

        alpha = np.min((np.max(((sd[i]-si3[i])/(sd[i]-si0[i]), 0)), 1))

        a[0, 0] = si0[i]-su[0, i]
        a[0, 4] = Aui[0]*(si0[i]-su[0, i])
        a[1, 1] = si0[i]-su[1, i]
        a[1, 4] = Aui[1]*(si0[i]-su[1, i])
        a[2, 0] = su[0, i]-si0[i]
        a[2, 1] = su[1, i]-si0[i]
        a[2, 2] = si0[i]-su[2, i]
        a[2, 4] = Aui[2]*(si0[i]-su[2, i])
        a[3, 2] = su[2, i]-si3[i]
        a[3, 3] = si3[i]-su[3, i]
        a[3, 5] = Aui[3]*(si3[i]-su[3, i])
        a[4, 2] = alpha*si3[i]+(1-alpha)*sd[i]-si0[i]
        a[4, 4] = Aui[0]*(su[0, i]-si0[i])+Aui[1]*(su[1, i]-si0[i])+Aui[1]*(su[1, i]-si0[i])
        a[4, 6] = Aid0*(sd[i]-si0[i])
        a[5, 3] = sp[i]-si3[i]
        a[5, 5] = Aui[3]*(su[3, i]-si3[i])
        a[6, 2] = (1-alpha)*(si3[i]-sd[i])
        a[6, 6] = Aid0*(si0[i]-sd[i])

        b = np.array([Vu[0]*dsudt[0, i]+F[0, i]*si0[i],
                      Vu[1]*dsudt[1, i]+F[1, i]*si0[i],
                      Vu[2]*dsudt[2, i]+F[2, i]*si0[i],
                      Vu[3]*dsudt[3, i]+F[3, i]*si0[i],
                      Vi0*dsi0dt[i]+(F[0, i]+F[1, i]+F[2, i])*(alpha*si3[i]+(1-alpha)*sd[i]-si0[i]),
                      Vi3*dsi3dt[i]+(F[0, i]+F[1, i]+F[2, i]+F[3, i])*(sp[i]-si3[i]),
                      Vd0*dsddt[i]+(F[0, i]+F[1, i]+F[2, i])*(1-alpha)*(si3[i]-sd[i]),
                     ])

        if heatflux == 1:
            # Temperature constrain
            a[7, 2] = alpha*ti3[i]+(1-alpha)*td[i]-ti0[i]
            a[7, 4] = Aui[0]*(tu[0, i]-ti0[i])+Aui[1]*(tu[1, i]-ti0[i])+Aui[1]*(tu[1, i]-ti0[i])
            a[7, 6] = Aid0*(td[i]-ti0[i])
            a[8, 3] = tp[i]-ti3[i]
            a[8, 5] = Aui[3]*(tu[3, i]-ti3[i])
            a[9, 2] = (1-alpha)*(ti3[i]-td[i])
            a[9, 6] = Aid0*(ti0[i]-td[i])

            # Temperature constrain
            bt = np.array([Vi0*dti0dt[i]+(F[0, i]+F[1, i]+F[2, i])*(alpha*ti3[i]+(1-alpha)*td[i]-ti0[i]),
                           Vi3*dti3dt[i]+(F[0, i]+F[1, i]+F[2, i]+F[3, i])*(tp[i]-ti3[i]),
                           Vd0*dtddt[i]+(F[0, i]+F[1, i]+F[2, i])*(1-alpha)*(ti3[i]-td[i]),
                          ])

            b = np.concatenate([b, bt])

        X[:, i] = lsq_linear(a, b, bounds=(np.zeros(7), np.ones(7)*np.inf)).x

    # ----------------------------------------------------------------------------------------------
    # Save data
    # svdata = 1
    if svpth != -1:
        import netCDF4 as nc
        # Write data into netCDF file
        print 'Saving data as netCDF4 file to '+svpth+'...'
        f = nc.Dataset(svpth+'inverse2.nc', 'w', format='NETCDF4')
        f.description = 'Outout of the inverse box model (boxMethod=1)'

        f.createDimension('time', None)
        f.createDimension('layer', 3)
        f.createDimension('box', 4)
        
        t_nc = f.createVariable('t', 'f8', ('time'))
        q02_nc = f.createVariable('Q02', 'f8', ('time'))
        q12_nc = f.createVariable('Q12', 'f8', ('time'))
        q23_nc = f.createVariable('Q23', 'f8', ('time'))
        q3p_nc = f.createVariable('Q3p', 'f8', ('time'))
        mix_ui0_nc = f.createVariable('Wui0', 'f8', ('time'))
        mix_ui3_nc = f.createVariable('Wui3', 'f8', ('time'))
        mix_id_nc = f.createVariable('Wid', 'f8', ('time'))
        s_nc = f.createVariable('s_box', 'f8', ('layer', 'box', 'time'))
        f_nc = f.createVariable('f_box', 'f8', ('box', 'time'))
        sp_nc = f.createVariable('sp', 'f8', ('time'))

        t_nc[:] = mt
        s_nc[:, :, :] = s
        f_nc[:, :] = F
        sp_nc[:] = sp
        q02_nc[:] = X[0, :]
        q12_nc[:] = X[1, :]
        q23_nc[:] = X[2, :]
        q3p_nc[:] = X[3, :]
        mix_ui0_nc[:] = X[4, :]
        mix_ui3_nc[:] = X[5, :]
        mix_id_nc[:] = X[6, :]

        f.close()

    inv = {'mt':    mt,
           's_box': s,
           'f_box': F,
           'si0':   si0,
           'sd':    sd,
           'sp':    sp,
           'Q02':   X[0, :],
           'Q12':   X[1, :],
           'Q23':   X[2, :],
           'Q3p':   X[3, :],
           'Wui0':  X[4, :],
           'Wui3':  X[5, :],
           'Wid':   X[6, :]
          }

    return inv


def box_inverse3(cal_clim=-1, svpth=-1, filterwindow=90):

    import numpy as np
    from box_gb.box_params import Aui, Aid, Vu, Vi, Vd, deltaSui, deltaSid, damp_mixing
    from box_gb.box_params import get_pacific
    from scipy.optimize import lsq_linear  # Solve linear system with constrains

    tidemixing = 1

    # Calculate/Read salinities and freshwater discharge in each box
    # The function box_avgbox_clim works fine on pool-97, but does not work on my laptop.
    # Thus use cal_clim = - 1 unless you want to change the averaging algorithm.
    if cal_clim == 1:
        from box_gb.box_ctd import box_avgbox_clim
        mt, s, F = box_avgbox_clim(svpth='./data/', filterwindow=filterwindow)
        mt, t, F = box_avgbox_clim(var='t', svpth='./data/', filterwindow=filterwindow)
    else:
        import netCDF4 as nc
        fh = nc.Dataset('./data/s_box.nc', mode='r')
        mt = fh.variables['t'][:]
        s = fh.variables['s_box'][:]
        F = fh.variables['f_box'][:]
        fh.close()
        
        fh = nc.Dataset('./data/t_box.nc', mode='r')
        t = fh.variables['t_box'][:]
        fh.close()

    su = s[0, :, :]
    si0 = (s[1, 0, :]*Vi[0]+s[1, 1, :]*Vi[1]+s[1, 2, :]*Vi[2])/(Vi[0]+Vi[1]+Vi[2])
    si3 = s[1, 3, :]
    sd = (s[2, 0, :]*Vd[0]+s[2, 1, :]*Vd[1]+s[2, 2, :]*Vd[2])/(Vd[0]+Vd[1]+Vd[2])
    Vi0 = Vi[0]+Vi[1]+Vi[2]
    Vi3 = Vi[3]
    Vd0 = Vd[0]+Vd[1]+Vd[2]
    Aid0 = Aid[0]+Aid[1]+Aid[2]
    # ----------------------------------------------------------------------------------------------
    dmt = np.diff(mt)*24*60*60  # [s]
    dsdt = np.diff(s, axis=2)/dmt  # [PSU s^-1]
    dsudt = dsdt[0, :, :]
    dsi0dt = (dsdt[1, 0, :]*Vi[0]+dsdt[1, 1, :]*Vi[1]+dsdt[1, 2, :]*Vi[2])/(Vi[0]+Vi[1]+Vi[2])
    dsi3dt = dsdt[1, 3, :]
    dsddt = (dsdt[2, 0, :]*Vd[0]+dsdt[2, 1, :]*Vd[1]+dsdt[2, 2, :]*Vd[2])/(Vd[0]+Vd[1]+Vd[2])
    # ----------------------------------------------------------------------------------------------
    mt = 0.5*(mt[:-1]+mt[1:])
    su = 0.5*(su[:, :-1]+su[:, 1:])  # [PSU]
    si0 = 0.5*(si0[:-1]+si0[1:])  # [PSU]
    si3 = 0.5*(si3[:-1]+si3[1:])  # [PSU]
    sd = 0.5*(sd[:-1]+sd[1:])  # [PSU]
    s = 0.5*(s[:, :, :-1]+s[:, :, 1:])  # [PSU]
    F = 0.5*(F[:, :-1]+F[:, 1:])  # [m^3s^-1]

    # Pacific Salinity (sp)
    # sp = get_pacific(mt, lat_sp=56.875, lon_sp=-137.625, filterwindow=filterwindow)
    sp = get_pacific(mt, filterwindow=filterwindow)  # [PSU]

    # ----------------------------------------------------------------------------------------------
    A = np.zeros((7*mt.size, 4+4*mt.size))
    B = np.zeros(7*mt.size)*np.nan

    for i in range(mt.size):
        # Solve linear system aX=b with constraint X>0
        a1 = np.zeros((7, 4))
        a2 = np.zeros((7, 4))
        c = np.zeros((4, 1))

        # --------------------------------------------------------------------------------
        # This version is created on 2016/09/01 -- see my notes
        # 2016/09/01
        # Combine intermediate and deep layer for box 0, 1, and 2

        alpha = np.min((np.max(((sd[i]-si3[i])/(sd[i]-si0[i]), 0)), 1))

        a1[0, 0] = si0[i]-su[0, i]
        a1[1, 1] = si0[i]-su[1, i]
        a1[2, 0] = su[0, i]-si0[i]
        a1[2, 1] = su[1, i]-si0[i]
        a1[2, 2] = si0[i]-su[2, i]
        a1[3, 2] = su[2, i]-si3[i]
        a1[3, 3] = si3[i]-su[3, i]
        a1[4, 2] = alpha*si3[i]+(1-alpha)*sd[i]-si0[i]
        a1[5, 3] = sp[i]-si3[i]
        a1[6, 2] = (1-alpha)*(si3[i]-sd[i])

        a2[0, 0] = Aui[0]*(si0[i]-su[0, i])
        a2[1, 0] = Aui[1]*(si0[i]-su[1, i])
        a2[2, 0] = Aui[2]*(si0[i]-su[2, i])
        a2[3, 1] = Aui[3]*(si3[i]-su[3, i])
        a2[4, 0] = Aui[0]*(su[0, i]-si0[i])+Aui[1]*(su[1, i]-si0[i])+Aui[1]*(su[1, i]-si0[i])
        a2[4, 2] = Aid0*(sd[i]-si0[i])
        a2[5, 1] = Aui[3]*(su[3, i]-si3[i])
        a2[5, 3] = sd[i]-si3[i]
        a2[6, 2] = Aid0*(si0[i]-sd[i])
        a2[6, 3] = si3[i]-sd[i]

        c[0] = su[2, i]-su[0, i]
        c[1] = su[2, i]-su[1, i]
        c[2] = su[3, i]-su[2, i]
        c[3] = su[3, i]-su[2, i]

        b = np.array([Vu[0]*dsudt[0, i]+F[0, i]*si0[i],
                      Vu[1]*dsudt[1, i]+F[1, i]*si0[i],
                      Vu[2]*dsudt[2, i]+F[2, i]*si0[i],
                      Vu[3]*dsudt[3, i]+F[3, i]*si0[i],
                      Vi0*dsi0dt[i]+(F[0, i]+F[1, i]+F[2, i])*(alpha*si3[i]+(1-alpha)*sd[i]-si0[i]),
                      Vi3*dsi3dt[i]+(F[0, i]+F[1, i]+F[2, i]+F[3, i])*(sp[i]-si3[i]),
                      Vd0*dsddt[i]+(F[0, i]+F[1, i]+F[2, i])*(1-alpha)*(si3[i]-sd[i]),
                     ])

        A[i*7:(i+1)*7, 0:4] = np.dot(a1, c)
        A[i*7:(i+1)*7, (i+1)*4:(i+2)*4] = a2
        B[i*7:(i+1)*7] = b

    X = lsq_linear(A, B, bounds=(np.zeros(4+4*mt.size), np.ones(4+4*mt.size)*np.inf)).x

    c0, c1, c2, c3 = X[0:4]
    w0 = X[4::4]
    w3 = X[5::4]
    wd = X[6::4]
    wk = X[7::4]

    # ----------------------------------------------------------------------------------------------
    # Save data
    # svdata = 1
    if svpth != -1:
        import netCDF4 as nc
        # Write data into netCDF file
        print 'Saving data as netCDF4 file to '+svpth+'...'
        f = nc.Dataset(svpth+'inverse3.nc', 'w', format='NETCDF4')
        f.description = 'Outout of the inverse box model (boxMethod=1)'

        f.createDimension('time', None)
        f.createDimension('layer', 3)
        f.createDimension('box', 4)
        f.createDimension('const', 1)
        
        t_nc = f.createVariable('t', 'f8', ('time'))
        c0_nc = f.createVariable('c0', 'f8', ('const'))
        c1_nc = f.createVariable('c1', 'f8', ('const'))
        c2_nc = f.createVariable('c2', 'f8', ('const'))
        c3_nc = f.createVariable('c3', 'f8', ('const'))
        mix_ui0_nc = f.createVariable('Wui0', 'f8', ('time'))
        mix_ui3_nc = f.createVariable('Wui3', 'f8', ('time'))
        mix_id_nc = f.createVariable('Wid', 'f8', ('time'))
        mix_hr_nc = f.createVariable('W3d', 'f8', ('time'))
        s_nc = f.createVariable('s_box', 'f8', ('layer', 'box', 'time'))
        f_nc = f.createVariable('f_box', 'f8', ('box', 'time'))
        sp_nc = f.createVariable('sp', 'f8', ('time'))

        t_nc[:] = mt
        s_nc[:, :, :] = s
        f_nc[:, :] = F
        sp_nc[:] = sp
        c0_nc[:] = c0
        c1_nc[:] = c1
        c2_nc[:] = c2
        c3_nc[:] = c3
        mix_ui0_nc[:] = w0
        mix_ui3_nc[:] = w3
        mix_id_nc[:] = wd
        mix_hr_nc[:] = wk

        f.close()

    inv = {'mt':    mt,
           's_box': s,
           'f_box': F,
           'si0':   si0,
           'sd':    sd,
           'sp':    sp,
           'c0':    c0,
           'c1':    c1,
           'c2':    c2,
           'c3':    c3,
           'Wui0':  w0,
           'Wui3':  w3,
           'Wid':   wd,
           'W3k':   wk
          }

    return inv
