import numpy as np
import hummingbirdParam as P


class ctrlPID:
    def __init__(self):
        # pitch

        tr_pitch = 1.0
        zeta_pitch = 0.707

        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        
        wn_pitch = 2.2 / tr_pitch
        self.kp_pitch = wn_pitch**2 / b_theta
        self.kd_pitch = 2 * zeta_pitch * wn_pitch / b_theta

        p_roots = np.roots([1, (2*zeta_pitch*wn_pitch) / b_theta, (wn_pitch**2) / b_theta])
        print('Pitch desired poles: ', p_roots)

        self.ki_pitch = 0.4

        self.integrator_theta = 0.0
        
        print('kp_pitch: ', self.kp_pitch)
        print('kd_pitch: ', self.kd_pitch) 
        print('ki_pitch: ', self.ki_pitch)
                
        # delayed variables
        self.theta_d1 = 0.
        self.theta_dot = 0.

        # roll
        # tuning parameters
        tr_roll = 0.3
        zeta_roll = 0.707
        
        # gain calculation
        b_phi = 1/P.J1x 
        
        #print('b_theta: ', b_theta)
        wn_roll = 2.2 / tr_roll
        self.kp_roll = wn_roll**2 / b_phi
        self.kd_roll = 2 * zeta_roll * wn_roll / b_phi
        
        # print gains to terminal
        print('kp_roll: ', self.kp_roll)
        #print('ki_roll: ', self.ki_roll)
        print('kd_roll: ', self.kd_roll) 
        
        # delayed variables
        self.phi_d1 = 0.
        self.phi_dot = 0.


        # yaw
        tr_yaw = 10 * tr_roll
        zeta_yaw = 0.707

        Fe = (P.g * (P.ell1 * P.m1 + P.ell2 * P.m2)) / P.ellT
        b_psi = (Fe * P.ellT) / (P.J1z + P.J2z + P.J3z + 
                                 (P.ell1**2 * P.m1) + 
                                 (P.ell2**2 * P.m2) + 
                                 (P.ell3x**2 * P.m3) + 
                                 (P.ell3y**2 * P.m3))        
        wn_yaw = 2.2 / tr_yaw
        self.kp_yaw = wn_yaw**2 / b_psi
        self.kd_yaw = 2 * zeta_yaw * wn_yaw / b_psi

        y_roots = np.roots([1, (2*zeta_yaw*wn_yaw) / b_psi, (wn_yaw**2)/ b_psi])
        print('Yaw desired poles: ', y_roots)

        self.ki_yaw = 0.001

        self.integrator_psi = 0.0

        # print gains to terminal
        print('kp_yaw: ', self.kp_yaw)
        print('kd_yaw: ', self.kd_yaw)
        print('ki_yaw: ', self.ki_yaw)

        self.psi_d1 = 0.
        self.psi_dot = 0.

        # sample rate of the controller
        self.Ts = P.Ts
        
        # dirty derivative parameters
        self.sigma = 0.05  # cutoff freq for dirty derivative
        self.beta = (2 * self.sigma - self.Ts) / (2 * self.sigma + self.Ts)
        

    def update(self, r: np.ndarray, y: np.ndarray):
        #############################################################################################
        # theta - longitudinal
        theta_ref = r[0, 0]
        theta = y[1, 0]
        force_fl = (P.m1 * P.ell1 + P.m2 * P.ell2) * P.g / P.ellT * np.cos(theta)  # feedforward force
        
        # compute errors
        error_theta = theta_ref - theta
        
        # update differentiators
        self.theta_dot = (self.beta * self.theta_dot + 
                          ((2.0 / (2.0 * self.sigma + self.Ts)) * (theta - self.theta_d1)))
        
        # update integrator with anti-windup
        if np.abs(self.theta_dot) <= 0.05:
            self.integrator_theta = self.integrator_theta + (P.Ts / 2) * (error_theta + (theta - self.theta_d1))

        # pitch control
        force = error_theta * self.kp_pitch - self.theta_dot * self.kd_pitch + self.ki_pitch * self.integrator_theta + force_fl
        force = force[0]

        # update all delayed variables
        self.theta_d1 = theta

        #############################################################################################

        # psi - lateral
        psi_ref = r[1, 0]
        psi = y[2, 0]

        error_psi = psi_ref - psi

        self.psi_dot = (self.beta * self.psi_dot +
                        ((2.0 / (2.0 * self.sigma + self.Ts)) * (psi - self.psi_d1)))

        # phi - lateral
        phi_ref = self.kp_yaw * error_psi - self.kd_yaw * self.psi_dot # comes from outer loop controller
        phi = y[0, 0]
        self.psi_d1 = psi

        
        # compute errors
        error_phi = phi_ref - phi
        
        # update differentiators
        self.phi_dot = (self.beta * self.phi_dot +
                        ((2.0 / (2.0 * self.sigma + self.Ts)) * (phi - self.phi_d1)))
        
        # update integrator with anti-windup
        if np.abs(self.psi_dot) <= 0.1:
            self.integrator_psi = self.integrator_psi + (P.Ts / 2) * (error_psi + (psi - self.psi_d1))

        # roll control - tau
        torque = error_phi * self.kp_roll - self.phi_dot * self.kd_roll + self.ki_yaw * self.integrator_psi
        torque = torque[0]

        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)

        # update all delayed variables
        self.phi_d1 = phi
        
        # return pwm plus reference signals
        phi_ref = phi_ref[0]
        return pwm, np.array([[phi_ref],[theta_ref], [psi_ref]])


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i, 0] > up_limit:
                u[i, 0] = up_limit
            if u[i, 0] < low_limit:
                u[i, 0] = low_limit
    return u