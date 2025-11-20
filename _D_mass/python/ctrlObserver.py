import numpy as np
import massParam as P
import control as cnt


class ctrlObserver:
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        zeta = 0.707
        omega_n = 2.2 / tr

        tr_obs = tr / 10.0
        zeta_obs = 0.707

        self.Pi = -1  # integrator pole

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

        des_char_poly = np.convolve([1, 2.0 * zeta * omega_n, omega_n**2], [1, -self.Pi])
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

        print("Poles: ",des_poles) 

        # observer design
        wn_obs = 2.2 / tr_obs
        des_obsv_char_poly = [1, 2 * zeta_obs * wn_obs, wn_obs**2]
        des_obsv_poles = np.roots(des_obsv_char_poly)

        if np.linalg.matrix_rank(cnt.ctrb(A.T, C.T)) != A.shape[0]:
            print("The system is not observable")
        else:
            self.L = cnt.place(A.T, C.T, des_obsv_poles).T
            print("L^T: ", self.L.T)

        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample
        self.x_hat = np.array([[0.0],   # z_hat
                               [0.0]])  # zdot_hat
        self.F_prev = 0.0  # previous force

    def update(self, z_r, y):
        x_hat = self.update_observer(y)
        # z = x[0, 0]
        error_z = z_r - x_hat[0, 0]

        # # integrate error with anti-windup
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        # # compute the state feedback controller
        F_unsat = -self.K @ x_hat - self.Ki * self.integrator

        # # compute total force
        F = saturate(F_unsat[0], P.F_max)
        self.F_prev = F

        return F, x_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat = self.x_hat + P.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)

        return self.x_hat
    
    def observer_f(self, x_hat, y_m):
        # xhatdot = A*xhat + B*F + L*(y - C*xhat)
        xhat_dot = P.Amat @ x_hat + P.Bmat * self.F_prev + self.L @ (y_m - P.Cmat @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u