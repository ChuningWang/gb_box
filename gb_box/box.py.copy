import numpy as np
from scipy.interpolate import interp1d
import pdb
from box_params import boxMethod, tide_sn, damp_mixing, damp_wind_mixing

if boxMethod==1:
    from box_params import t, f, ts, Sp, lambda_u, lambda_i, lambda_d, Rui, Rid, R12, R23, R3p, Rt, deltaSui, deltaSid
elif boxMethod==2:
    from box_params import t, f, ts, Sp, lambda_u, lambda_i, lambda_d, Rui, Rid, R23u, R23i, R12, R2p, R3p, Rt

if damp_wind_mixing==1:
    from box_params import awnd, awnd_crit

# ------------------------------------------------------------
if boxMethod==1:

    def box_ode(mt, y):

        ft = interp1d(t,f.T)(mt)

        S = np.reshape(y, (4, 3)).T

        if damp_mixing == 1:
            # Damp mixing if stratification is too strong
            Rui_d = np.zeros(Rui.shape)
            Rid_d = np.zeros(Rid.shape)
            for i in range(Rui.size):
                if S[1, i] >= S[0, i]:
                    Rui_d[i] = Rui[i]*np.min([deltaSui[i]/(S[1, i]-S[0, i]), 1])
                else:
                    Rui_d[i] = Rui[i]
                if S[2, i] >= S[1, i]:
                    Rid_d[i] = Rid[i]*np.min([deltaSid[i]/(S[2, i]-S[1, i]), 1])
                else:
                    Rid_d[i] = Rid[i]
        else:
            Rui_d, Rid_d = Rui, Rid

        if damp_wind_mixing==1:
            # Damp surface-intermediate mixing when wind speed is low
            # This is the beta test
            w = interp1d(t, awnd)(mt)
            if w<awnd_crit:
                Rui_d[3] = 0

        if tide_sn==1:
            # Add a spring-neap tidal cycle
            Rui_d[3] = Rui_d[3]*(1+0.4*(np.sin(2*np.pi*mt/(365/24))-1))

        q02 = S[0, 2]-S[0, 0]
        q02i = q02-ft[0]

        q12 = R12*(S[0, 2]-S[0, 1])
        q12i = q12-ft[1]

        q23 = R23*(S[0, 3]-S[0, 2])
        q2iu = q23-q02-q12-ft[2]
        q23i = q02i+q12i+q2iu

        q3p = R3p*(Sp-S[0, 3])
        # q3p = R3p*q23+ft[3]
        q3iu = q3p-q23-ft[3]
        q3pi = q23i+q3iu

        dS = np.zeros((3, 4))

        dS[0, 0] = (1./lambda_u[0]/ts)*(
                                        +q02i*S[1, 0]-q02*S[0, 0]
                                        +Rui_d[0]*(S[1, 0]-S[0, 0])
                                       )
        
        dS[1, 0] = (1./lambda_i[0]/ts)*(
                                        +q02i*(S[1, 2]-S[1, 0])
                                        -Rui_d[0]*(S[1, 0]-S[0, 0])+Rid_d[0]*(S[2, 0]-S[1, 0])
                                       )
        
        dS[2, 0] = (1./lambda_d[0]/ts)*(
                                        -Rid_d[0]*(S[2, 0]-S[1, 0])
                                       )

        dS[0, 1] = (1./lambda_u[1]/ts)*(
                                        +q12i*S[1, 1]-q12*S[0, 1]
                                        +Rui_d[1]*(S[1, 1]-S[0, 1])
                                       )

        dS[1, 1] = (1./lambda_i[1]/ts)*(
                                        +q12i*(S[1, 2]-S[1, 1])
                                        -Rui_d[1]*(S[1, 1]-S[0, 1])+Rid_d[1]*(S[2, 1]-S[1, 1])
                                       )

        dS[2, 1] = (1./lambda_d[1]/ts)*(
                                        -Rid_d[1]*(S[2, 1]-S[1, 1])
                                       )

        dS[0, 2] = (1./lambda_u[2]/ts)*(
                                        +q02*S[0, 0]+q12*S[0, 1]+q2iu*S[1, 2]-q23*S[0, 2]
                                        +Rui_d[2]*(S[1, 2]-S[0, 2])
                                       )

        dS[1, 2] = (1./lambda_i[2]/ts)*(
                                        +q23i*(S[1, 3]-S[1, 2])
                                        -Rui_d[2]*(S[1, 2]-S[0, 2])+Rid_d[2]*(S[2, 2]-S[1, 2])
                                       )

        dS[2, 2] = (1./lambda_d[2]/ts)*(
                                        -Rid_d[2]*(S[2, 2]-S[1, 2])
                                       )

        dS[0, 3] = (1./lambda_u[3]/ts)*(
                                        +q23*S[0, 2]+q3iu*S[1, 3]-q3p*S[0, 3]
                                        +Rui_d[3]*(S[1, 3]-S[0, 3])
                                       )

        dS[1, 3] = (1./lambda_i[3]/ts)*(
                                        +q3pi*(Sp-S[1, 3])
                                        -Rui_d[3]*(S[1, 3]-S[0, 3])+Rid_d[3]*(S[2, 3]-S[1, 3])
                                       )+Rt*(Sp-S[1, 3])

        dS[2, 3] = (1./lambda_d[3]/ts)*(
                                        -Rid_d[3]*(S[2, 3]-S[1, 3])
                                       )

        dy = np.reshape(dS.T, 12)
        # pdb.set_trace()
        return dy

elif boxMethod==2:

    def box_ode(mt, y):

        # idx = np.where(t==mt)[0][0]
        # f0 = f[idx, 0]
        # f1 = f[idx, 1]
        # f2 = f[idx, 2]
        # f3 = f[idx, 3]
        f_t = interp1d(t,f.T)(mt)

        S = np.reshape(y, (4, 3)).T

        q03 = S[0, 3]-S[0, 0]
        q12 = S[0, 2]-S[0, 1]
        q3p = Sp-S[0, 3]
        q2p = Sp-S[0, 2]

        dS = np.zeros((3, 4))

        dS[0, 0] = (1./lambda_u[0]/ts)*(
                                        +(q03-f_t[0])*S[1, 0]-q03*S[0, 0]
                                        +Rui[0]*(S[1, 0]-S[0, 0])
                                       )
        
        dS[1, 0] = (1./lambda_i[0]/ts)*(
                                        +(q03-f_t[0])*(S[1, 3]-S[1, 0])
                                        -Rui[0]*(S[1, 0]-S[0, 0])+Rid[0]*(S[2, 0]-S[1, 0])
                                       )
        
        dS[2, 0] = (1./lambda_d[0]/ts)*(
                                        -Rid[0]*(S[2, 0]-S[1, 0])
                                       )

        dS[0, 1] = (1./lambda_u[1]/ts)*(
                                        +(R12*q12-f_t[1])*S[1, 1]-R12*q12*S[0, 1]
                                        +Rui[1]*(S[1, 1]-S[0, 1])
                                       )

        dS[1, 1] = (1./lambda_i[1]/ts)*(
                                        +(R12*q12-f_t[1])*(S[1, 2]-S[1, 1])
                                        -Rui[1]*(S[1, 1]-S[0, 1])+Rid[1]*(S[2, 1]-S[1, 1])
                                       )

        dS[2, 1] = (1./lambda_d[1]/ts)*(
                                        -Rid[1]*(S[2, 1]-S[1, 1])
                                       )

        dS[0, 2] = (1./lambda_u[2]/ts)*(
                                        +R12*q12*S[0, 1]+(R2p*q2p-R12*q12-f_t[2])*S[1, 2]-R2p*q2p*S[0, 2]
                                        +Rui[2]*(S[1, 2]-S[0, 2])
                                        +R23u*(S[0, 3]-S[0, 2])
                                       )

        dS[1, 2] = (1./lambda_i[2]/ts)*(
                                        +(R2p*q2p-f_t[1]-f_t[2])*(Sp-S[1, 2])
                                        -Rui[2]*(S[1, 2]-S[0, 2])+Rid[2]*(S[2, 2]-S[1, 2])
                                        +R23i*(S[1, 3]-S[1, 2])
                                       )+Rt*(Sp-S[1, 2])

        dS[2, 2] = (1./lambda_d[2]/ts)*(
                                        -Rid[2]*(S[2, 2]-S[1, 2])
                                       )

        dS[0, 3] = (1./lambda_u[3]/ts)*(
                                        +q03*S[0, 0]+(R3p*q3p-q03-f_t[3])*S[1, 3]-R3p*q3p*S[0, 3]
                                        +Rui[3]*(S[1, 3]-S[0, 3])
                                        -R23u*(S[0, 3]-S[0, 2])
                                       )

        dS[1, 3] = (1./lambda_i[3]/ts)*(
                                        +(R3p*q3p-f_t[0]-f_t[3])*(Sp-S[1, 3])
                                        -Rui[3]*(S[1, 3]-S[0, 3])+Rid[3]*(S[2, 3]-S[1, 3])
                                        -R23i*(S[1, 3]-S[1, 2])
                                       )+Rt*(Sp-S[1, 3])

        dS[2, 3] = (1./lambda_d[3]/ts)*(
                                        -Rid[3]*(S[2, 3]-S[1, 3])
                                       )

        dy = np.reshape(dS.T, 12)
        # pdb.set_trace()
        return dy

# ------------------------------------------------------------

def rk_solver(y0, mt, b, output=1):

    mt = mt[:]

    dy = b(mt[0], y0)
    y = np.zeros((dy.size, mt.size))
    y[:, 0] = y0

    ct = 0
    for i in range(1, mt.size):
        ct = ct+1
        h = mt[i]-mt[i-1]
        k1 = b(mt[i-1],     y[:, i-1])
        k2 = b(mt[i-1]+h/2, y[:, i-1]+k1*h/2)
        k3 = b(mt[i-1]+h/2, y[:, i-1]+k2*h/2)
        k4 = b(mt[i-1]+h,   y[:, i-1]+k3*h)
        y[:, i] = y[:, i-1]+(1./6)*h*(k1+2*k2+2*k3+k4)
        # pdb.set_trace()
        if output==1:
            if np.remainder(ct, 10)==0:
                print 'S='
                print y[:, i]

    return y

# ------------------------------------------------------------

