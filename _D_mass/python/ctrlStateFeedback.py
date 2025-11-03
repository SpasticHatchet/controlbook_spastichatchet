import numpy as np
import massParam as P


class ctrlPID:
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        zeta = 0.707

        # compute PD gains
        # open loop char polynomial and poles
        a1 = P.b / P.m
        a0 = P.k / P.m
        wn = 2.2/tr
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2

        # dirty derivative gains
        self.sigma = 0.05  
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)  # dirty derivative gain

        #----------------------------------------------------------
        # variables for integrator and differentiator
        self.z_dot = P.zdot0             # estimated derivative of y
        self.z_prev = P.z0              # Signal y delayed by one sample
        self.error_dot = 0.0         # estimated derivative of error
        self.error_prev = 0.0          # Error delayed by one sample
        # self.integrator = 0.0        # integrator
                
    def update(self, z_r, y):
        z = y[0][0]
        # Compute the current error
        error = z_r - z

        # integrator anti - windup
        if self.z_dot < 0.1:
            # integrate error
            self.integrator = self.integrator + (P.Ts / 2) * (error + self.error_prev)

        # differentiate z
        self.z_dot = (2 * self.sigma - P.Ts) / (2 * self.sigma + P.Ts) * self.z_dot + \
                            (2.0 / (2.0*self.sigma + P.Ts)) * (z - self.z_prev)
        
        # PID control
        F_unsat = self.kp * error + self.ki * self.integrator - self.kd * self.z_dot
        F = saturate(F_unsat, P.F_max)
        
        # update delayed variables
        self.error_prev = error
        self.z_prev = z
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u