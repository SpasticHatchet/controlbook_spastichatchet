import numpy as np
import massParam as P
#from massDynamics import saturate


class ctrlPD:
    def __init__(self):
        #  tuning parameters
        tr = 2 # tuned for faster rise time before saturation.
        zeta = 0.7

        # desired natural frequency
        wn = 2.2 / tr
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2

        a1 = P.b / P.m
        a0 = P.k / P.m
        b0 = 1 / P.m

        self.Kd = (alpha1 - a1) / b0
        self.Kp = (alpha0 - a0) / b0

        # self.Kd = 4.5
        # self.Kp = 12.0

        print('kp: ', self.Kp)
        print('kd: ', self.Kd)

    def update(self, r, state):

        z = state[0,0]
        zdot = state[1,0]

        F = (self.Kp * (r - z)) - (self.Kd * zdot)

        if np.abs(F) > 6:
            print('F :', F)

        #F = saturate(F, P.F_max)
        return F