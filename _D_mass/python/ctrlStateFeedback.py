import numpy as np
import massParam as P
import control as cnt


class ctrlStateFeedback:
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        zeta = 0.707
        omega_n = 2.2 / tr

        A = P.Amat
        B = P.Bmat
        Cr = P.Cmat[0]

        des_char_poly = [1, 2.0 * zeta * omega_n, omega_n**2]
        des_poles = np.roots(des_char_poly)

        n = np.size(A, 1)

        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != n:
            print("The system is not controllable")
        else:
            self.K = cnt.place(A, B, des_poles)
            print(self.K)
            A_BK = np.linalg.inv(A - B @ self.K)
            self.Kr = -1.0 / (Cr @ A_BK @ B)

        print('K: ', self.K)
        print('Kr: ', self.Kr)
        print(des_poles)

        # # compute PD gains
        # # open loop char polynomial and poles
        # a1 = P.b / P.m
        # a0 = P.k / P.m
        # wn = 2.2/tr
        # alpha1 = 2.0 * zeta * wn
        # alpha0 = wn**2

        # # dirty derivative gains
        # self.sigma = 0.05  
        # self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)  # dirty derivative gain

        # #----------------------------------------------------------
        # # variables for integrator and differentiator
        # self.z_dot = P.zdot0             # estimated derivative of y
        # self.z_prev = P.z0              # Signal y delayed by one sample
        # self.error_dot = 0.0         # estimated derivative of error
        # self.error_prev = 0.0          # Error delayed by one sample
        # # self.integrator = 0.0        # integrator
                
    # def update(self, z_r, y):
    #     z = y[0][0]
    #     # Compute the current error
    #     error = z_r - z

    #     # integrator anti - windup
    #     if self.z_dot < 0.1:
    #         # integrate error
    #         self.integrator = self.integrator + (P.Ts / 2) * (error + self.error_prev)

    #     # differentiate z
    #     self.z_dot = (2 * self.sigma - P.Ts) / (2 * self.sigma + P.Ts) * self.z_dot + \
    #                         (2.0 / (2.0*self.sigma + P.Ts)) * (z - self.z_prev)
        
    #     # PID control
    #     F_unsat = self.kp * error + self.ki * self.integrator - self.kd * self.z_dot
    #     F = saturate(F_unsat, P.F_max)
        
    #     # update delayed variables
    #     self.error_prev = error
    #     self.z_prev = z
    #     return F

    # def update(self, theta_r, x):
    #     theta = x[0, 0]
    #     # compute feedback linearizing torque tau_fl
    #     tau_fl = P.m * P.g * (P.ell / 2.0) * np.cos(theta)

    #     # Compute the state feedback controller
    #     tau_tilde = -self.K @ x + self.kr * theta_r

    #     # compute total torque
    #     tau = saturate(tau_fl + tau_tilde[0, 0], P.tau_max)

    #     return tau

    def update(self, z_r, x):
        # compute the state feedback controller
        F_tilde = -self.K @ x + self.Kr * z_r

        # compute total force
        F = saturate(F_tilde[0, 0], P.F_max)

        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u