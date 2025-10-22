# This is modeled from Design Study A, ctrlPDhw7.py

import numpy as np
import massParam as P
from massDynamics import saturate

class ctrlPD:
    def __init__(self):
        # Declare pole locations
        #des_poles_CE = np.convolve([1, 3], [1, 4]) # desired poles from (s+3)(s+4)
        p1 = -1
        p2 = -1.5
        a1 = P.b / P.m
        a0 = P.k / P.m
        b0 = 1 / P.m
        # alpha0 =
        # alpha1 =
        self.Kp = 1.5*P.m - P.k
        self.Kd = 2.5*P.m - P.b

        # PD gains
        print('Kp: ', self.Kp)
        print('Kd: ', self.Kd)

        #self.use_feedback_linearization = True

    def update(self, r, state):

        z = state[0,0]
        zdot = state[1,0]

        F = (self.Kp * (r - z)) - (self.Kd * zdot)
        #u = (P.Kp * (r - z)) - (P.Kd * zdot)
        F = saturate(F, P.F_max)
        return F
    
    #def update(self, theta_r, x):
        # theta = x[0, 0]
        # thetadot = x[1, 0]

        # # compute the linearized torque using PD control
        # tau_tilde = self.kp * (theta_r - theta) - self.kd * thetadot

        # # Include feadback linearization torque or equilibrium 
        # # torque around theta_e
        # if self.use_feedback_linearization:
        #     # Feedback linearized torque
        #     tau_fl = P.m * P.g * (P.ell / 2.0) * np.cos(theta)

        #     # compute total torque
        #     tau = tau_fl + tau_tilde
        # else:
        #     # equilibrium torque around theta_e = 0
        #     theta_e = 0.0
        #     tau_e = P.m * P.g * P.ell / 2.0 * np.cos(theta_e)

        #     # compute total torque
        #     tau = tau_e + tau_tilde

        # return tau

