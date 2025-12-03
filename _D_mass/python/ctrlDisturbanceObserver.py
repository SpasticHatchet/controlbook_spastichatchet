import numpy as np
import massParam as P
import control as cnt


class ctrlDisturbanceObserver:
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        zeta = 0.707
        omega_n = 2.2 / tr

        tr_obs = tr / 10.0
        zeta_obs = 0.707

        integrator_pole = -1.0  # integrator pole
        disturbance_pole = -6.0  # disturbance observer pole

        A = P.Amat
        B = P.Bmat
        C = P.Cmat
        
        # build augmented A1 = [[A, 0],
        #                       [-Cr, 0]]
        Cr_row = np.reshape(np.asarray(C[0]), (1, A.shape[1]))  # ensure shape (1,2)
        A1 = np.vstack((
            np.hstack((A, np.zeros((A.shape[0], 1)))),
            np.hstack((-Cr_row, np.zeros((1, 1))))
        ))
        B1 = np.vstack((B, np.zeros((1, 1))))

        print("A1:", A1)
        print("B1:", B1)

        des_char_poly = np.convolve([1, 2.0 * zeta * omega_n, omega_n**2], 
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)

        n = np.size(A1, 1)

        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != n:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0, 0:2]
            self.Ki = K1[0, 2]
            print("K: ",self.K)
            print("Ki: ",self.Ki)
            # A_BK = np.linalg.inv(A1 - B1 @ K1)
            # self.Kr = -1.0 / (Cr @ A_BK @ B1)

        #print("Poles: ",des_poles) 

        # disturbance observer design
        A2 = np.block([[A, B],
                       [np.zeros((1, A.shape[1])), 0]])
        print("A2:", A2)
        B2 = B1 # same as B1
        C2 = np.block([C, np.zeros((C.shape[0], 1))])
        print("C2:", C2)

        des_obsv_char_poly = np.convolve([1, 2 * zeta_obs * omega_n, omega_n**2],
                                         [1, -disturbance_pole])
        des_obsv_poles = np.roots(des_obsv_char_poly)
        # Compute the gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != A2.shape[0]:
            print("The system is not observerable")
        else:
            L2 = cnt.place(A2.T, C2.T, des_obsv_poles).T
        print('L^T: ', L2.T)

        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample
        self.obsv_state = np.array([
            [0.0],  # z_hat_0
            [0.0],  # zdot_hat_0
            [0.0],  # estimate of the disturbance
        ])
        self.force_d1 = 0.0  # previous force
        self.L = L2
        self.A = A2
        self.B = B1
        self.C = C2

    def update(self, z_r, y_m):
        x_hat, d_hat = self.update_observer(y_m)
        z_hat = x_hat[0, 0]
        error_z = z_r - z_hat

        # integrate error with anti-windup
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        # compute the state feedback controller
        force_unsat = -self.K @ x_hat - self.Ki * self.integrator - d_hat

        # compute total force
        force = saturate(force_unsat[0], P.F_max)
        self.force_d1 = force

        return force, x_hat, d_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.obsv_state, y_m)
        F2 = self.observer_f(self.obsv_state + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.obsv_state + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.obsv_state + P.Ts * F3, y_m)
        self.obsv_state = self.obsv_state + P.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)
        x_hat = self.obsv_state[0:2]
        d_hat = self.obsv_state[2, 0]

        return x_hat, d_hat

    def observer_f(self, x_hat, y_m):
        xhat_dot = self.A @ x_hat + self.B * self.force_d1 + self.L @ (y_m - self.C @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u