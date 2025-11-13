import numpy as np
import control as cnt
import blockbeamParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------

        A = P.Amat
        B = P.Bmat
        Cr = P.Cmat[0]

        # print("A:", A)
        # print("B:", B)
        # print("Cr:", Cr)
        
        # build augmented matrices
        Cr_row = np.reshape(np.asarray(Cr), (1, A.shape[1]))  # ensure shape
        # print("Cr_row:", Cr_row)

        A1 = np.vstack((
            np.hstack((A, np.zeros((A.shape[0], 1)))),
            np.hstack((-Cr_row, np.zeros((1, 1))))
        ))
        # print("A1:", A1)

        B1 = np.vstack((B, np.zeros((1, 1))))
        # print("B1:", B1)
        
        tr_theta = 0.25
        tr_z = 10 * tr_theta
        omegan_z = 2.2 / tr_z # 2.2/tr
        zeta_z = 0.707
        omegan_theta = 2.2 / tr_theta
        zeta_theta = 0.707

        self.Pi = -12  # integrator pole

        des_char_poly = np.convolve( np.convolve(
            [1, 2 * zeta_z * omegan_z, omegan_z**2],
            [1, 2 * zeta_theta * omegan_theta, omegan_theta**2]),
            [1, -self.Pi])


        des_poles = np.roots(des_char_poly)

        C_AB = cnt.ctrb(A1, B1)

        n = np.size(A1, 1)
        rank = np.linalg.matrix_rank(C_AB)
        detC_AB = np.linalg.det(C_AB)
        #%% 
        # print(n)
        # print(rank)
        # print(detC_AB)

        K1 = cnt.place(A1, B1, des_poles)
        self.K = K1[0, 0:4]
        self.Ki = K1[0, 4]
        print("K: ",self.K)
        print("Ki: ",self.Ki)

        print("Poles: ",des_poles)

        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample   
        # %%

    
    def update(self, z_r, x):
        z = x[0,0]
        error_z = z_r - z

        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        F_unsat = -self.K @ x - self.Ki * self.integrator

        # x_tilde = x - np.array([[P.ze], [0], [0], [0]])
        # zr_tilde = z_r - P.ze   # P.ze is the same as C_r*x_e 

        # # equilibrium force
        # F_e = P.m1*P.g*(P.ze/P.length) + P.m2*P.g/2.0

        # # Compute the state feedback controller
        # F_tilde = -self.K @ x_tilde + self.Ki * zr_tilde
        # F_unsat = F_e + F_tilde
        F = saturate(F_unsat[0], P.F_max)

        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u