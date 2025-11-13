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

        C_AB = cnt.ctrb(A, B)

        n = np.size(A, 1)
        rank = np.linalg.matrix_rank(C_AB)
        detC_AB = np.linalg.det(C_AB)
        #%% 
        # print(n)
        # print(rank)
        # print(detC_AB)

        self.K = cnt.place(A, B, des_poles)
        print("K:", self.K)

        A_BK = A - B @ self.K
        self.Kr = -1.0 / (C @ np.linalg.inv(A_BK) @ B)

        print("Kr:", self.Kr)
        # %%

    
    def update(self, z_r, x):
        z = x[0,0]
        x_tilde = x - np.array([[P.ze], [0], [0], [0]])
        zr_tilde = z_r - P.ze   # P.ze is the same as C_r*x_e 

        # equilibrium force
        F_e = P.m1*P.g*(P.ze/P.length) + P.m2*P.g/2.0

        # Compute the state feedback controller
        F_tilde = -self.K @ x_tilde + self.Kr * zr_tilde
        F_unsat = F_e + F_tilde
        F = saturate(F_unsat[0,0], P.F_max)

        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
