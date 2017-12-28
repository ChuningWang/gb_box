""" ODE function for the box model. """

import numpy as np
import pdb

class BoxSolver(object):
    """ solver for the box model. """

    def __init__(self):

        self.solution = {}
        self.deltaS = 1e-4
        self.max_loop = 100

        return None

    def __call__(self, box):

        box_num = len(box.boxes_list)*3
        # generate initial condition
        sinit = box.consts['sinit']*np.ones(box_num)

        if box.info['clim']:
            print 'Climatology run...'
            s_end = np.zeros(box_num)
            cts = 0
            while np.max(np.abs(s_end - sinit)) > self.deltaS:
                if cts > 0:
                    sinit = self.solution[:, -1]

                self.box_solver(box, sinit)
                s_end = self.solution[:, -1]
                cts += 1
                print 'Loop ' + str(cts)
                if cts >= self.max_loop:
                    print 'Maximum loop reached, terminating'
                    break

        return None

    def box_solver(self, box, sinit):
        """ solve the box model. """

        # parse constants, model parameters.
        class KeyVars(object):
            pass
        KeyVars.consts = box.consts
        KeyVars.volumes = box.volumes
        KeyVars.areas = box.areas
        KeyVars.rivers = box.BoxRivers.data
        KeyVars.sp = box.BoxSp.data

        # solve the ODEs
        if box.info['box_method'] == 1:
            self.solution = rk_solver(sinit, box.time, box_ode, KeyVars, verb=False)
        elif box.info['box_method'] == 2:
            self.solution = rk_solver(sinit, box.time, box_ode2, KeyVars, verb=False)

        return None


def box_ode(time, y, kv):
    """ ODE functions for box method 1. """

    s = np.reshape(y, (3, 3)).T

    ft0 = np.interp(time, kv.rivers['time'], kv.rivers['river0'])
    ft1 = np.interp(time, kv.rivers['time'], kv.rivers['river1'])
    ft2 = np.interp(time, kv.rivers['time'], kv.rivers['river2'])
    sp = np.interp(time, kv.sp['time'], kv.sp['salt'])

    # damp mixing/transport when stratification is high
    # add a spring-neap tidal cycle
    omega0ui_damp = \
        abs(min([kv.consts['deltasui0']/kv.consts['s0']/(s[1, 0]-s[0, 0]), 1])) * \
        0.5 * (1 + np.sin(2 * np.pi * time / (365 / 24)))
    omega0id_damp = \
        0.5*(1+np.sign(s[2, 0]-s[1, 0]))*1 + 0.5*(1-np.sign(s[2, 0]-s[1, 0]))*1.e3
    omega1ui_damp = \
        abs(min([kv.consts['deltas1']/kv.consts['s0']/(s[1, 1]-s[0, 1]), 1])) * \
        0.5 * (1 + np.sin(2 * np.pi * time / (365 / 24)))

    c0_damp = \
        min([kv.consts['deltasc0']/kv.consts['s0']/abs(s[1, 0]-s[0, 0]), 1]) * \
        (1 + 0.25*np.sin(2 * np.pi * time / (365 / 24)))

    # apply damping to c0
    c0brs = kv.consts['c0brs']*c0_damp
    c1brs = kv.consts['c1brs']

    # convert discharge from m^3s^-1 to m^3day^-1
    ft0 = ft0*kv.consts['s2d']
    ft1 = ft1*kv.consts['s2d']
    ft2 = ft2*kv.consts['s2d']

    ds = np.zeros(y.shape)
    q0 = s[0, 1] - s[0, 0]
    q0i = q0 - ft0/c0brs
    q1 = s[0, 2] - s[0, 1]
    q1i = q1 - ft0/c1brs \
             - ft1/c1brs
    q2 = q1 + ft2/c1brs

    ds[0] = 1./kv.volumes['Vu0']*(
        c0brs*(-q0*s[0, 0] + q0i*s[1, 0])
        + omega0ui_damp*kv.consts['omegasi0']*kv.areas['Asi0']*(s[1, 0] - s[0, 0])
                                 )

    ds[1] = 1./kv.volumes['Vi0']*(
        c0brs*(-q0i*s[1, 0] + q0i*s[1, 1])
        - omega0ui_damp*kv.consts['omegasi0']*kv.areas['Asi0']*(s[1, 0] - s[0, 0])
        + omega0id_damp*kv.consts['omegaid0']*kv.areas['Aid0']*(s[2, 0] - s[1, 0])
                                 )

    ds[2] = 1./kv.volumes['Vd0']*(
        - omega0id_damp*kv.consts['omegaid0']*kv.areas['Aid0']*(s[2, 0] - s[1, 0])
                                 )

    ds[3] = 1./kv.volumes['Vu1']*(
        c1brs*(-q1*s[0, 1] + q1i*s[1, 1])
        + kv.consts['c0brs']*(-q0i*s[0, 1] + q0*s[0, 0])
        + omega1ui_damp*kv.consts['omegasi1']*kv.areas['Asi1']*(s[1, 1] - s[0, 1])
                                 )

    ds[4] = 1./kv.volumes['Vi1']*(
        c1brs*(-q1i*s[1, 1] + q1i*s[1, 2])
        + kv.consts['c0brs']*(-q0i*s[1, 1] + q0i*s[0, 1])
        - omega1ui_damp*kv.consts['omegasi1']*kv.areas['Asi1']*(s[1, 1] - s[0, 1])
                                 )

    ds[6] = 1./kv.volumes['Vu2']*(
        c1brs*(-q2*s[0, 2] + q1*s[0, 1])
        + kv.consts['omegasi2']*kv.areas['Asi2']*(s[1, 2] - s[0, 2])
                                 )

    ds[7] = 1./kv.volumes['Vi2']*(
        c1brs*(-q1i*s[1, 2] + q1i*sp/kv.consts['s0'])
        - kv.consts['omegasi2']*kv.areas['Asi2']*(s[1, 2] - s[0, 2])
                                 ) \
        + 1./kv.consts['tr']*(sp/kv.consts['s0'] - s[1, 2])

    return ds

def box_ode2(time, y, kv):
    """
    ODE functions for box method 2.

    2017-12-27
    Add a mixing term between box 3 and box 4.
    """

    # ------------------- preparations ----------------------------------------
    # assign values
    su = y[0::3]
    si = y[1::3]
    sd = y[2::3]

    # interpolate river discharge
    ft = np.zeros(5)
    for i in range(5):
        ft[i] = np.interp(time, kv.rivers['time'], kv.rivers['river' + str(i)])
        # convert from m^3s^-1 to m^3day^-1
        ft[i] = ft[i]*kv.consts['s2d']

    # interpolate pacific salinity
    sp = np.interp(time, kv.sp['time'], kv.sp['salt'])

    # ------------------- add-on features -------------------------------------
    b1 = 0.5
    q_sn = 1. + b1*np.sin(2*np.pi*time/(365/24))
    q_damp = np.zeros(5)
    for i in range(5):
        q_damp[i] = kv.consts['deltasq']/(kv.consts['s0']*max([si[i]-su[i], 1.e-5]))
        q_damp[i] = min([q_damp[i], 1])
    # mixing damping/amplification
    mui = np.zeros(5)
    mid = np.zeros(5)
    for i in range(5):
        mui[i] = kv.consts['deltasui']/(kv.consts['s0']*max([si[i]-su[i], 1.e-7]))
        mui[i] = min([max([mui[i], 1]), 100])
        mid[i] = kv.consts['deltasid']/(kv.consts['s0']*max([sd[2]-si[i], 1.e-7]))
        mid[i] = min([max([mid[i], 1]), 100])

    # ------------------- calculate volume transports -------------------------
    q0 = q_sn*q_damp[0]*kv.consts['c0brs']*(su[2] - su[0])
    q0e = q0 - ft[0]
    q1 = q_sn*q_damp[1]*kv.consts['c1brs']*(su[2] - su[1])
    q1e = q1 - ft[1]
    q2 = q_sn*q_damp[2]*kv.consts['c2brs']*(su[3] - su[2])
    q2e = q2 - q0 - q1 - ft[2]
    q2i = q0e + q1e + q2e
    q3 = q2*kv.consts['c_amp'] + ft[3]
    q3e = q3 - q2 - ft[3]
    q3i = q3e + q2i
    q4 = q3 + ft[4]
    q4i = q3i

    ds = np.zeros(y.shape)

    ds[0] = 1./kv.volumes['Vu0']*(
        (-q0*su[0] + q0e*si[0])
        + kv.consts['omegav']*kv.areas['Aui0']*mui[0]*(si[0] - su[0])
                                 )

    ds[1] = 1./kv.volumes['Vi0']*(
        (-q0e*si[0] + q0e*si[2])
        - kv.consts['omegav']*kv.areas['Aui0']*mui[0]*(si[0] - su[0])
        + kv.consts['omegav']*kv.areas['Aid0']*mid[0]*(sd[2] - si[0])
                                 )

    ds[3] = 1./kv.volumes['Vu1']*(
        (-q1*su[1] + q1e*si[1])
        + kv.consts['omegav']*kv.areas['Aui1']*mui[1]*(si[1] - su[1])
                                 )

    ds[4] = 1./kv.volumes['Vi1']*(
        (-q1e*si[1] + q1e*si[2])
        - kv.consts['omegav']*kv.areas['Aui1']*mui[1]*(si[1] - su[1])
        + kv.consts['omegav']*kv.areas['Aid1']*mid[1]*(sd[2] - si[1])
                                 )

    ds[6] = 1./kv.volumes['Vu2']*(
        (-q2*su[2] + q0*su[0] + q1*su[1] + q2e*si[2])
        + kv.consts['omegav']*kv.areas['Aui2']*mui[2]*(si[2] - su[2])
                                 )

    ds[7] = 1./kv.volumes['Vi2']*(
        (-q2i*si[2] + q2i*si[3])
        - kv.consts['omegav']*kv.areas['Aui2']*mui[2]*(si[2] - su[2])
        + kv.consts['omegav']*kv.areas['Aid2']*mid[2]*(sd[2] - si[2])
                                 )

    ds[8] = 1./kv.volumes['Vd2']*(
        - kv.consts['omegav']*kv.areas['Aid0']*mid[0]*(sd[2] - si[0])
        - kv.consts['omegav']*kv.areas['Aid1']*mid[1]*(sd[2] - si[1])
        - kv.consts['omegav']*kv.areas['Aid2']*mid[2]*(sd[2] - si[2])
                                 )

    ds[9] = 1./kv.volumes['Vu3']*(
        (-q3*su[3] + q2*su[2] + q3e*si[3])
        + kv.consts['omegav2']*kv.areas['Aui3']*mui[3]*(si[3] - su[3])
        + kv.consts['omegah']*(su[4] - su[3])
                                 )

    ds[10] = 1./kv.volumes['Vi3']*(
        (-q2i*si[3] - q3e*si[3] + q3i*si[4])
        - kv.consts['omegav2']*kv.areas['Aui3']*mui[3]*(si[3] - su[3])
        + kv.consts['omegah']*(si[4] - si[3])
                                  )

    ds[12] = 1./kv.volumes['Vu4']*(
        (-q4*su[4] + q3*su[3])
        + kv.consts['omegav2']*kv.areas['Aui4']*mui[4]*(si[4] - su[4])
        - kv.consts['omegah']*(su[4] - su[3])
                                  )

    ds[13] = 1./kv.volumes['Vi4']*(
        (-q3i*si[4] + q4i*sp/kv.consts['s0'])
        - kv.consts['omegav2']*kv.areas['Aui4']*mui[4]*(si[4] - su[4])
        - kv.consts['omegah']*(si[4] - si[3])
                                  ) \
        + 1./kv.consts['tr']*(sp/kv.consts['s0'] - si[4])

    return ds

def rk_solver(y0, time, odefunc, kvars={}, verb=True):
    """ solve ODE using 4th order Runge-Kutta method. """

    dy = odefunc(time[0], y0, kvars)
    y = np.zeros((len(dy), len(time)))
    y[:, 0] = y0

    ct = 0
    for i in range(1, len(time)):
        ct = ct+1
        h = time[i] - time[i-1]
        k1 = odefunc(time[i-1],     y[:, i-1], kvars)
        k2 = odefunc(time[i-1]+h/2, y[:, i-1]+k1*h/2, kvars)
        k3 = odefunc(time[i-1]+h/2, y[:, i-1]+k2*h/2, kvars)
        k4 = odefunc(time[i-1]+h,   y[:, i-1]+k3*h, kvars)
        y[:, i] = y[:, i-1] + (1./6)*h*(k1+2*k2+2*k3+k4)

        if verb:
            if np.remainder(ct, 10) == 0:
                print 'y='
                print y[:, i]

    return y
