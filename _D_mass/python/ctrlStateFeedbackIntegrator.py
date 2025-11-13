import numpy as np
import massParam as P
import control as cnt


class ctrlStateFeedbackIntegrator:
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        zeta = 0.707
        omega_n = 2.2 / tr

        A = P.Amat
        B = P.Bmat
        Cr = P.Cmat[0]

        # build augmented A1 = [[A, 0],
        #                       [-Cr, 0]]
        Cr_row = np.reshape(np.asarray(Cr), (1, A.shape[1]))  # ensure shape (1,2)
        A1 = np.vstack((
            np.hstack((A, np.zeros((A.shape[0], 1)))),
            np.hstack((-Cr_row, np.zeros((1, 1))))
        ))
        B1 = np.vstack((B, np.zeros((1, 1))))

        print("A1:", A1)
        print("B1:", B1)

        self.Pi = -1  # integrator pole

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
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample        

    def update(self, z_r, x):
        z = x[0, 0]
        error_z = z_r - z

        # integrate error with anti-windup
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)

        
        self.error_d1 = error_z
        # compute the state feedback controller
        F_unsat = -self.K @ x - self.Ki * self.integrator

        # compute total force
        F = saturate(F_unsat[0], P.F_max)

        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u