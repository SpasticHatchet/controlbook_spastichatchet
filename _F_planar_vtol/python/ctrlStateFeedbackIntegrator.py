
import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedback:
    def __init__(self):

        self.Pi_lon = -1.0
        self.Pi_lat = -1.0

        self.F_max = P.F_max

        # self.Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force 

        tr_z = 2.0
        omegan_z = 2.2 / tr_z
        zeta_z = 0.707

        tr_theta = tr_z / 10.0
        omegan_th = 2.2 / tr_theta
        zeta_th = 0.707

        tr_h = 2.0
        omegan_h = 2.2 / tr_h
        zeta_h = 0.707

        A_lon = P.Amat_lon
        B_lon = P.Bmat_lon
        C_lon = P.Cmat_lon[0]
        D_lon = P.Dmat_lon
        A_lat = P.Amat_lat
        B_lat = P.Bmat_lat
        C_lat = P.Cmat_lat[0]
        D_lat = P.Dmat_lat

        A1_lon = np.vstack((np.hstack((A_lon, np.zeros((A_lon.shape[0],1)))),
                            np.hstack((-C_lon, np.array([[0.0]])))))
        B1_lon = np.vstack((B_lon, np.array([[0.0]])))
        C1_lon = np.hstack((C_lon, np.array([[0.0]])))

        A1_lat = np.vstack((np.hstack((A_lat, np.zeros((A_lat.shape[0],1)))),
                            np.hstack((-C_lat, np.array([[0.0]])))))
        B1_lat = np.vstack((B_lat, np.array([[0.0]])))
        C1_lat = np.hstack((C_lat, np.array([[0.0]])))

        # gain calculation
        des_char_poly_lon = np.convolve([1.0, 2.0 * zeta_h * omegan_h, omegan_h ** 2],
                                        [1.0, -self.Pi_lon])
        des_poles_lon = np.roots(des_char_poly_lon)
        des_char_poly_lat = np.convolve( np.convolve([1.0, 2.0 * zeta_z * omegan_z, omegan_z ** 2],
                                        [1.0, 2.0 * zeta_th * omegan_th, omegan_th ** 2]),
                                        [1.0, -self.Pi_lat])
        des_poles_lat = np.roots(des_char_poly_lat)
        
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A_lon, B_lon)) != A_lon.shape[0]:
            print("The longitudinal system is not controllable")
        else:
            K1_lon = cnt.place(A_lon, B_lon, des_poles_lon)
            self.K_lon = K1_lon[0:2]
            self.Ki_lon = K1_lon[0, 2]
        #     self.kr_lon = -1.0 / (C_lon @ np.linalg.inv(A_lon - B_lon @ self.K_lon) @ B_lon)
        
        if np.linalg.matrix_rank(cnt.ctrb(A_lat, B_lat)) != A_lat.shape[0]:
            print("The lateral system is not controllable")
        else:
            K1_lat = cnt.place(A_lat, B_lat, des_poles_lat)
            self.K_lat = K1_lat[0:4]
            self.Ki_lat = K1_lat[0, 4]
        #     self.kr_lat = -1.0 / (Cr_lat @ np.linalg.inv(A_lat - B_lat @ self.K_lat) @ B_lat)

        print('K_lon: ', self.K_lon)
        print('Ki_lon: ', self.Ki_lon)
        print('K_lat: ', self.K_lat)
        print('Ki_lat: ', self.Ki_lat)

    def update(self, r, x):
        z_r = r[0,0]
        h_r = r[1,0]
        z = x[0,0]
        h = x[1,0]
        theta = x[2,0]

        # Construct the states
        x_lon = np.array([[x[1,0]], [x[4,0]]])
        x_lat = np.array([[x[0,0]], [x[2,0]], [x[3,0]], [x[5,0]]])

        # Compute the state feedback controllers (F and tau) 
        F_tilde = -self.K_lon @ x_lon + self.kr_lon * h_r
        F = self.Fe + F_tilde[0,0]
        tau_tilde = -self.K_lat @ x_lat + self.kr_lat * z_r
        
        # Saturating the motor thrusts for left and right motors
        motor_thrusts = P.mixing @ np.array([[F], [tau_tilde[0,0]]])
        motor_thrusts[0,0] = saturate(motor_thrusts[0,0], P.F_max)
        motor_thrusts[1,0] = saturate(motor_thrusts[1,0], P.F_max)

        return motor_thrusts

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

# import numpy as np
# import VTOLParam as P
# import control as cnt

# class ctrlStateFeedback:
#     def __init__(self):
#         # tuning parameters for each set of poles
#         # for longitudinal dynamics
#         tr_h = 2.0
#         wn_h = 2.2 / tr_h
#         zeta_h = 0.707
        
#         # for lateral dynamics
#         tr_z = 2.0
#         wn_z = 2.2 / tr_z
#         zeta_z = 0.707
#         M = 10.0
#         tr_th = tr_z / M
#         wn_th = 2.2 / tr_th
#         zeta_th = 0.707

#         # It is important to note that we formulate the dynamics as longitudinal dynamics
#         # and lateral dynamics, but this is not strictly necessary. It may help us to see
#         # and understand the effect of the two parts that are mostly decoupled (for the
#         # linearized system at least). However, it is totally acceptable to also form
#         # one large set of state space equations that includes A_lon and A_lat in one
#         # A matrix (similary for B, C, etc.). See example commented out below in this 
#         # "init" function.

#         # State Space Equations for lateral and longitudinal dynamics
#         self.Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force 

#         # for longitudinal dynamics
#         A_lon = np.array([[0.0, 1.0],
#                         [0.0, 0.0]])
#         B_lon = np.array([[0.0],
#                         [1.0 / (P.mc + 2.0 * P.mr)]])
#         C_lon = np.array([[1.0, 0.0]])

#         # for lateral dynamics
#         A_lat = np.array([[0.0, 0.0, 1.0, 0.0],
#                         [0.0, 0.0, 0.0, 1.0],
#                         [0.0, -self.Fe / (P.mc + 2.0 * P.mr), -(P.mu / (P.mc + 2.0 * P.mr)), 0.0],
#                         [0.0, 0.0, 0.0, 0.0]])
#         B_lat = np.array([[0.0],
#                         [0.0],
#                         [0.0],
#                         [1.0 / (P.Jc + 2 * P.mr * P.d ** 2)]])
#         C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
#                         [0.0, 1.0, 0.0, 0.0]])
        
#         # gain calculation
#         des_char_poly_lon = [1.0, 2.0 * zeta_h * wn_h, wn_h ** 2]
#         des_poles_lon = np.roots(des_char_poly_lon)
#         des_char_poly_lat = np.convolve([1.0, 2.0 * zeta_z * wn_z, wn_z ** 2],
#                                         [1.0, 2.0 * zeta_th * wn_th, wn_th ** 2])
#         des_poles_lat = np.roots(des_char_poly_lat)
        
#         # Compute the gains if the system is controllable
#         if np.linalg.matrix_rank(cnt.ctrb(A_lon, B_lon)) != A_lon.shape[0]:
#             print("The longitudinal system is not controllable")
#         else:
#             self.K_lon = cnt.place(A_lon, B_lon, des_poles_lon)
#             self.kr_lon = -1.0 / (C_lon @ np.linalg.inv(A_lon - B_lon @ self.K_lon) @ B_lon)
        
#         if np.linalg.matrix_rank(cnt.ctrb(A_lat, B_lat)) != A_lat.shape[0]:
#             print("The lateral system is not controllable")
#         else:
#             self.K_lat = cnt.place(A_lat, B_lat, des_poles_lat)
#             Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])
#             self.kr_lat = -1.0 / (Cr_lat @ np.linalg.inv(A_lat - B_lat @ self.K_lat) @ B_lat)

#         print('K_lon: ', self.K_lon)
#         print('kr_lon: ', self.kr_lon)
#         print('K_lat: ', self.K_lat)
#         print('kr_lat: ', self.kr_lat)

#         # this is how we could do it if we wanted to use one full state space model:
#         # #state is defined in this order (same order as dynamics file)=> 
#         #       ['z', 'h', 'theta', 'z_dot', 'h_dot', 'theta_dot']

#         # des_polls_all = np.concatenate((des_poles_lon, des_poles_lat))

#         # # these models are hand-typed here, but they could come directly from SymPy as well. 
#         # A_full = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
#         #                    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
#         #                    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
#         #                    [0.0, 0.0, -self.Fe / (P.mc + 2.0 * P.mr), -P.mu / (P.mc + 2.0 * P.mr),  0.0, 0.0],
#         #                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#         #                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

#         # B_full = np.array([[0.0, 0.0],
#         #                    [0.0, 0.0],
#         #                    [0.0, 0.0],
#         #                    [0.0, 0.0],
#         #                    [1.0 / (P.mc + 2.0 * P.mr), 0.0],
#         #                    [0.0, 1.0 / (P.Jc + 2.0 * P.mr * P.d ** 2)]])

#         # if np.linalg.matrix_rank(cnt.ctrb(A_full, B_full)) != A_full.shape[0]:
#         #     print("The full system is not controllable")
#         # else:
#         #     C_full = np.array([[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
#         #                     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
#         #     self.K_full = cnt.place(A_full, B_full, des_polls_all)
#         #     self.kr_full = -np.linalg.inv(C_full @ np.linalg.inv(A_full - B_full @ self.K_full) @ B_full)

#         # print("Full K matrix: \n", self.K_full)
#         # print("Full kr matrix: \n", self.kr_full)

#     def update(self, r, x):
#         z_r = r[0,0]
#         h_r = r[1,0]
#         z = x[0,0]
#         h = x[1,0]
#         theta = x[2,0]

#         # Construct the states
#         x_lon = np.array([[x[1,0]], [x[4,0]]])
#         x_lat = np.array([[x[0,0]], [x[2,0]], [x[3,0]], [x[5,0]]])

#         # Compute the state feedback controllers (F and tau) 
#         F_tilde = -self.K_lon @ x_lon + self.kr_lon * h_r
#         F = self.Fe + F_tilde[0,0]
#         tau_tilde = -self.K_lat @ x_lat + self.kr_lat * z_r
        
#         # Saturating the motor thrusts for left and right motors
#         motor_thrusts = P.mixing @ np.array([[F], [tau_tilde[0,0]]])
#         motor_thrusts[0,0] = saturate(motor_thrusts[0,0], P.F_max)
#         motor_thrusts[1,0] = saturate(motor_thrusts[1,0], P.F_max)

#         return motor_thrusts


# def saturate(u, limit):
#     if abs(u) > limit:
#         u = limit*np.sign(u)
#     return u

