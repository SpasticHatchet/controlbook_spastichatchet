import numpy as np
import blockbeamParam as P

class ctrlPD:
    def __init__(self):
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        # tr_th = 0.5         # Ries time for part a) - inner loop
        tr_th = 1.0       # Rise time for inner loop (theta)
        zeta_th = 0.707       # inner loop Damping Coefficient

        # saturation limits
        #self.F_max = 5             		  # Max Force, N
        # error_max = 1        		  # Max step size,m
        # theta_max = 30.0 * np.pi / 180.0  # Max theta, rads
        A = ((P.m2 * P.length**2.0) / 3.0) + (P.m1 * (P.length / 2.0)**2.0)
        #---------------------------------------------------
        #                    Inner Loop
        #---------------------------------------------------
        # parameters of the open loop transfer function
        b0_th = P.length / A
        a0_th = 0.0
        a1_th = 0.0

        # coefficients for desired inner loop
        wn_th = 2.2 / tr_th     # Natural frequency

        # compute gains
        self.kp_th = (wn_th**2 + a0_th) / b0_th
        self.kd_th = (2.0 * zeta_th * wn_th) / b0_th

        DC_gain_th = b0_th * self.kp_th / (b0_th * self.kp_th + a0_th)

        #---------------------------------------------------
        #                    Outer Loop
        #---------------------------------------------------
        # coefficients for desired outer loop
        M = 10.0           # Time scale separation 
        zeta_z = 0.707     # outer loop Damping Coefficient
        tr_z = M * tr_th   # desired rise time, s
        wn_z = 2.2 / tr_z  # desired natural frequency

        b0_z = -P.g
        a0_z = 0.0
        a1_z = 0.0

        # compute gains
        #a = wn_z**2*np.sqrt(2.0*P.length/3.0/P.g)-2.0*zeta_z*wn_z
        self.kp_z = (wn_z**2 + a0_z) / b0_z
        self.kd_z = (2.0 * zeta_z * wn_z) / b0_z
        

        # print control gains to terminal        
        print('DC_gain', DC_gain_th)
        print('kp_th: ', self.kp_th)
        print('kd_th: ', self.kd_th)
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)

        #---------------------------------------------------
        #                    zero canceling filter
        #---------------------------------------------------
        # self.filter = zeroCancelingFilter(DC_gain_th)

    def update(self, z_r, state):
        z = state[0, 0]
        theta = state[1, 0]
        zdot = state[2, 0]
        thetadot = state[3, 0]

        # the reference angle for theta comes from the
        # outer loop PD control
        theta_r = self.kp_z * (z_r - z) - self.kd_z * zdot

        # the force applied to the cart comes from the
        # inner loop PD control
        F = self.kp_th * (theta_r - theta) - self.kd_th * thetadot
        Fe = (P.g / P.length) * (P.m2 * (P.length / 2.0) + P.m1*z) * np.cos(theta)
        F_total = F + Fe
        F_sat = F_total #saturate(F_total, P.F_max)

        return F_sat
    
def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u