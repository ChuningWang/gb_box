""" solve the ODEs using a Runge-Kutta solver. """

import numpy
from gb_box import boxPrep

class BoxSolver(object, box):

    def rk_solver(y0, time, odefunc, kv={}, print_verb=True):
        """ solve ODE using 4th order Runge-Kutta method. """

        dy = odefunc(time[0], y0, kv)
        y = np.zeros((len(dy), len(time)))
        y[:, 0] = y0

        ct = 0
        for i in range(1, len(time)):
            ct = ct+1
            h = time[i] - time[i-1]
            k1 = odefunc(time[i-1],     y[:, i-1], kv)
            k2 = odefunc(time[i-1]+h/2, y[:, i-1]+k1*h/2, kv)
            k3 = odefunc(time[i-1]+h/2, y[:, i-1]+k2*h/2, kv)
            k4 = odefunc(time[i-1]+h,   y[:, i-1]+k3*h, kv)
            y[:, i] = y[:, i-1] + (1./6)*h*(k1+2*k2+2*k3+k4)

            if print_verb:
                if np.remainder(ct, 10) == 0:
                    print 'S='
                    print y[:, i]

        return y
