import numpy as np
import control as cnt
import blockbeamParam as P

class ctrlStateFeedback:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------

        A = P.Amat
        B = P.Bmat
        C = P.Cmat[0]
        
        tr_theta = 0.25
        tr_z = 10 * tr_theta
        omegan_z = 2.2 / tr_z # 2.2/tr
        zeta_z = 0.707
        omegan_theta = 2.2 / tr_theta
        zeta_theta = 0.707

        des_char_poly = np.convolve(
            [1, 2 * zeta_z * omegan_z, omegan_z**2],
            [1, 2 * zeta_theta * omegan_theta, omegan_theta**2])


        des_poles = np.roots(des_char_poly)

        C_AB = cnt.ctrb(P.Amat, P.Bmat)

        n = np.size(P.Amat, 1)
        rank = np.linalg.matrix_rank(C_AB)
        detC_AB = np.linalg.det(C_AB)
        #%% 
        # print(n)
        # print(rank)
        # print(detC_AB)

        self.K = cnt.place(P.Amat, P.Bmat, des_poles)
        print("K:")
        print(self.K)

        A_BK = P.Amat - P.Bmat @ self.K
        self.Kr = -1.0 / (P.Cmat[0] @ np.linalg.inv(A_BK) @ P.Bmat)

        print("Kr:")
        print(self.Kr)
        # %%

    def update(self, z_r, x):
        z = x[0, 0]
        theta = x[1, 0]

        # compute feedback linearizing force F_fl
        # This cancels the gravity terms in the thetaddot equation
        F_fl = (P.m1 * P.g * z + P.m2 * P.g * P.length / 2.0) / P.length

        # Compute the state feedback controller
        # After feedback linearization, apply control directly to states
        F_tilde = -self.K @ x + self.Kr * z_r

        # compute total force
        F = saturate(F_fl + F_tilde[0, 0], P.F_max)

        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

