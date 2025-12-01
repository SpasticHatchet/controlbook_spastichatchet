import numpy as np
from scipy import signal
import control as cnt
import blockbeamParam as P


class ctrlObserver:
    def __init__(self):
        # tuning parameters
        tr_z = 2.0        # rise time for position
        tr_theta = 1.0    # rise time for angle
        zeta_z = 0.95     # damping ratio position
        zeta_th = 0.95    # damping ratio angle
        integrator_pole = -5.0

        # pick observer poles
        tr_z_obs = tr_z / 10.0
        tr_th_obs = tr_theta / 10.0
        
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        self.A = np.array([
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, -P.g, 0.0, 0.0],
            [-P.m1*P.g/((P.m2*P.length**2)/3.0+P.m1*(P.length/2.0)**2), \
             0.0, 0.0, 0.0]])
        self.B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length / (P.m2 * P.length ** 2 / 3.0 \
                            + P.m1 * P.length ** 2 / 4.0)]])
        self.C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0]])
        self.Cr = np.array([[1.0, 0.0, 0.0, 0.0]])  

        # form augmented system for integral control
        A1 = np.vstack((
                np.hstack((self.A, np.zeros((4,1)))),
                np.hstack((-self.Cr, np.zeros((1,1))))))
        print("A1:",A1)
        B1 = np.vstack((self.B, np.zeros((1,1))))
        
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

        # compute observer gains
        wn_z_obs = 0.5*np.pi/(tr_z_obs*np.sqrt(1-zeta_z**2))
        wn_th_obs = 0.5*np.pi/(tr_th_obs*np.sqrt(1-zeta_th**2))
        des_obs_char_poly = np.convolve(
            [1, 2*zeta_z*wn_z_obs, wn_z_obs**2],
            [1, 2*zeta_th*wn_th_obs, wn_th_obs**2])
        des_obs_poles = np.roots(des_obs_char_poly)

        # Compute the gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != self.A.T.shape[0]:
            print("The system is not observable")
        else:
            # .T transposes the output
            self.L = cnt.place(self.A.T, self.C.T, des_obs_poles).T
            print("A shape:",self.A.shape)
            print(("C shape,",self.C.shape))
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('L^T: ', self.L.T)

        # variables to implement integrator
        self.integrator_z = 0.0  # integrator
        self.error_z_prev = 0.0  # error signal delayed by 1 sample

        # initializing estimated state variables
        self.x_hat = np.array([
            [0.0],  # initial estimate for z_hat
            [0.0],  # initial estimate for theta_hat
            [0.0],  # initial estimate for z_hat_dot
            [0.0]])  # initial estimate for theta_hat_dot
        self.F_prev = 0.0  # Computed Force, delayed by one sample
        self.Ts = P.Ts

    def update(self, z_r, y):
        # update the observer and extract z_hat
        x_hat = self.update_observer(y)
        z_hat = self.Cr @ x_hat

        # integrate error
        error_z = z_r - z_hat
        self.integrator_z = self.integrator_z \
            + (P.Ts / 2.0) * (error_z + self.error_z_prev)
        self.error_z_prev = error_z

        # Construct the linearized state
        xe = np.array([[P.ze], [0.0], [0.0], [0.0]])
        x_tilde = x_hat - xe

        # equilibrium force
        F_e = P.m1*P.g*P.ze/P.length + P.m2*P.g/2.0

        # Compute the state feedback controller
        F_tilde = -self.K @ x_tilde \
                  - self.ki * self.integrator_z
        F_unsat = F_e + F_tilde
        F = saturate(F_unsat[0][0], P.F_max)
        self.integratorAntiWindup(F, F_unsat)
        self.F_prev = F

        return F, x_hat

    def update_observer(self, y):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y)
        F2 = self.observer_f(self.x_hat + self.Ts / 2 * F1, y)
        F3 = self.observer_f(self.x_hat + self.Ts / 2 * F2, y)
        F4 = self.observer_f(self.x_hat + self.Ts * F3, y)
        self.x_hat += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y):
        # xhatdot = A*(xhat-xe) + B*u + L(y-C*xhat)
        xe = np.array([[P.ze], [0.0], [0.0], [0.0]])

        # equilibrium force
        F_e = P.m1*P.g*P.ze/P.length + P.m2*P.g/2.0

        # observer dynamics
        xhat_dot = self.A @ (x_hat - xe) \
                   + self.B * (self.F_prev - F_e) \
                   + self.L @ (y - self.C @ x_hat)
        
        return xhat_dot

    def integratorAntiWindup(self, F, F_unsat):
        # integrator anti - windup
        if self.ki != 0.0:
            self.integrator_z = self.integrator_z \
                              + P.Ts/self.ki*(F-F_unsat)


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

