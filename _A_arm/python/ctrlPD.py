import numpy as np
import armParam as P


class ctrlPD:
    def __init__(self):
        #  tuning parameters
        tr = 0.37 # tuned for faster rise time before saturation.
        zeta = 0.7

        # desired natural frequency
        wn = 2.2 / tr
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2

        # compute PD gains
        self.kp = alpha0*(P.m * P.ell**2) / 3.0
        self.kd = (P.m * P.ell**2) \
                    / 3.0 * (alpha1 - 3.0 * P.b / (P.m * P.ell**2))
        print('kp: ', self.kp)
        print('kd: ', self.kd)

    def update(self, r, state):

            z = state[0,0]
            zdot = state[1,0]

            F = (self.Kp * (r - z)) - (self.Kd * zdot)
            #u = (P.Kp * (r - z)) - (P.Kd * zdot)
            F = saturate(F, P.F_max)
            return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








