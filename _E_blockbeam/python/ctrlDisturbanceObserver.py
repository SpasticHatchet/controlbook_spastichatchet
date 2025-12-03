import numpy as np
from scipy import signal
import control as cnt
import blockbeamParam as P

class ctrlDisturbanceObserver:
    def __init__(self):
        # tuning parameters
        tr_z = 2.0        # rise time for position
        tr_theta = 1.0    # rise time for angle
        zeta_z = 0.95     # damping ratio position
        zeta_th = 0.95    # damping ratio angle
        integrator_pole = -5.0
        dist_observer_pole = -10.0

        # pick observer poles
        tr_z_obs = tr_z / 10.0
        tr_th_obs = tr_theta / 10.0
        
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = P.Amat
        B = P.Bmat
        C = P.Cmat
        Cr = C[0]

        # form augmented system for integral control
        A1 = np.vstack((
                np.hstack((A, np.zeros((4,1)))),
                np.hstack((-Cr, np.zeros((1,1))))))
        print("A1:",A1)
        B1 = np.vstack((B, np.zeros((1,1))))
        
        # gain calculation
        wn_th = 0.5*np.pi/(tr_theta*np.sqrt(1-zeta_th**2))
        wn_z = 0.5*np.pi/(tr_z*np.sqrt(1-zeta_z**2))
        des_char_poly = np.convolve(
            np.convolve([1, 2*zeta_z*wn_z, wn_z**2],
                        [1, 2*zeta_th*wn_th, wn_th**2]),
            [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != A1.shape[0]:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]

        # disturbance observer design
        A2 = np.block([[A, B],
                       [np.zeros((1, A.shape[1])), 0]])
        print("A2:",A2)
        B2 = B1  # same as B1
        C2 = np.block([C, np.zeros((C.shape[0], 1))])
        print("C2:",C2)


        wn_z_obs = 0.5*np.pi/(tr_z_obs*np.sqrt(1-zeta_z**2))
        wn_th_obs = 0.5*np.pi/(tr_th_obs*np.sqrt(1-zeta_th**2))
        des_obs_char_poly = np.convolve(np.convolve(
            [1, 2*zeta_z*wn_z_obs, wn_z_obs**2],
            [1, 2*zeta_th*wn_th_obs, wn_th_obs**2]),
            [1, -dist_observer_pole])
        des_obs_poles = np.roots(des_obs_char_poly)

        # Compute the gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != A2.T.shape[0]:
            print("The system is not observable")
        else:
            # .T transposes the output
            L2 = cnt.place(A2.T, C2.T, des_obs_poles).T
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('L^T: ', L2.T)

        # variables to implement integrator
        self.integrator_z = 0.0  # integrator
        self.error_z_d1 = 0.0  # error signal delayed by 1 sample
        self.obsv_state = np.array([
            [0.0],  # z_hat_0
            [0.0],  # theta_hat_0
            [0.0],  # z_hat_dot_0
            [0.0],  # theta_hat_dot_0
            [0.0]])  # estimate of the disturbance
        self.force_d1 = 0.0  # Computed Force, delayed by one sample
        self.L = L2
        self.A = A2
        self.B = B1
        self.C = C2
        self.Cr = Cr
        self.Ts = P.Ts

    def update(self, z_r, y_m):
        x_hat, d_hat = self.update_observer(y_m)
        z_hat = self.Cr @ x_hat

        # integrate error
        error_z = z_r - z_hat
        self.integrator_z = self.integrator_z \
            + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z

        # Construct the linearized state
        xe = np.array([[P.ze], [0.0], [0.0], [0.0]])
        x_tilde = x_hat - xe

        # equilibrium force
        force_e = P.m1*P.g*P.ze/P.length + P.m2*P.g/2.0

        # Compute the state feedback controller
        force_tilde = -self.K @ x_tilde \
                  - self.ki * self.integrator_z \
                  - d_hat
        force_unsat = force_e + force_tilde
        force = saturate(force_unsat[0][0], P.F_max)
        self.integratorAntiWindup(force, force_unsat)
        self.force_d1 = force

        return force, x_hat, d_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.obsv_state, y_m)
        F2 = self.observer_f(self.obsv_state + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.obsv_state + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.obsv_state + self.Ts * F3, y_m)
        self.obsv_state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        x_hat = self.obsv_state[0:4]
        d_hat = self.obsv_state[4, 0]

        return x_hat, d_hat

    def observer_f(self, x_hat, y_m):
        # xhatdot = A*(xhat-xe) + B*u + L(y-C*xhat)
        xe = np.array([[P.ze], [0.0], [0.0], [0.0], [0.0]])

        # equilibrium force
        force_e = P.m1*P.g*P.ze/P.length + P.m2*P.g/2.0

        # observer dynamics
        xhat_dot = self.A @ (x_hat - xe) \
                   + self.B * (self.force_d1 - force_e) \
                   + self.L @ (y_m - self.C @ x_hat)
        
        return xhat_dot

    def integratorAntiWindup(self, force, force_unsat):
        # integrator anti - windup
        if self.ki != 0.0:
            self.integrator_z = self.integrator_z \
                              + P.Ts/self.ki*(force-force_unsat)


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

