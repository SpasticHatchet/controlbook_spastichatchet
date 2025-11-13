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
            # print(self.K)
            A_BK = np.linalg.inv(A - B @ self.K)
            self.Kr = -1.0 / (Cr @ A_BK @ B)

        print('K: ', self.K)
        print('Kr: ', self.Kr)
        print(des_poles)

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