
import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):

        self.Pi_lon = -1.0
        self.Pi_lat = -1.0

        self.F_max = P.F_max

        self.Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force 

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
        print("A1_lon:", A1_lon)
        B1_lon = np.vstack((B_lon, np.array([[0.0]])))
        C1_lon = np.hstack((C_lon, np.array([[0.0]])))

        A1_lat = np.vstack((np.hstack((A_lat, np.zeros((A_lat.shape[0],1)))),
                            np.hstack((-C_lat, np.array([[0.0]])))))
        print("A1_lat:", A1_lat)
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
        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != A1_lon.shape[0]:
            print("The longitudinal system is not controllable")
        else:
            K1_lon = cnt.place(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0, 0:2]  # Extract first 2 gains from row vector
            self.Ki_lon = K1_lon[0, 2]
            print("K1_lon:", K1_lon)
            print("K_lon:", self.K_lon)
            print("Ki_lon:", self.Ki_lon)
        
        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != A1_lat.shape[0]:
            print("The lateral system is not controllable")
        else:
            K1_lat = cnt.place(A1_lat, B1_lat, des_poles_lat)
            self.K_lat = K1_lat[0, 0:4]  # Extract first 4 gains from row vector
            self.Ki_lat = K1_lat[0, 4]
            print("K1_lat:", K1_lat)
            print("K_lat:", self.K_lat)
            print("Ki_lat:", self.Ki_lat)

        print('K_lon: ', self.K_lon)
        print('Ki_lon: ', self.Ki_lon)
        print('K_lat: ', self.K_lat)
        print('Ki_lat: ', self.Ki_lat)

        self.integrator_lon = 0.0
        self.integrator_lat = 0.0
        self.error_d1_lon = 0.0
        self.error_d1_lat = 0.0

    def update(self, r, x):
        z_r = r[0,0]
        h_r = r[1,0]

        z = x[0,0]
        h = x[1,0]

        error_z = z_r - z
        error_h = h_r - h

        # Longitudinal (h) integrator - integrates h error
        self.integrator_lon = self.integrator_lon + (P.Ts / 2.0) * (error_h + self.error_d1_lon)
        self.error_d1_lon = error_h

        # Lateral (z) integrator - integrates z error
        self.integrator_lat = self.integrator_lat + (P.Ts / 2.0) * (error_z + self.error_d1_lat)
        self.error_d1_lat = error_z

        # Construct the state vectors
        # Longitudinal: [h, hdot] from state vector [z, h, theta, zdot, hdot, thetadot]
        x_lon = np.array([[x[1,0]], [x[4,0]]])
        # Lateral: [z, theta, zdot, thetadot]
        x_lat = np.array([[x[0,0]], [x[2,0]], [x[3,0]], [x[5,0]]])

        # Longitudinal control (F) - uses longitudinal states and h integrator
        F_unsat = -self.K_lon @ x_lon - self.Ki_lon * self.integrator_lon + self.Fe
        F_unsat = F_unsat[0]
        # Lateral control (tau) - uses lateral states and z integrator
        tau_unsat = -self.K_lat @ x_lat - self.Ki_lat * self.integrator_lat
        tau_unsat = tau_unsat[0]

        motor_thrusts = P.mixing @ np.array([[F_unsat], [tau_unsat]])
        motor_thrusts[0,0] = saturate(motor_thrusts[0,0], self.F_max)
        motor_thrusts[1,0] = saturate(motor_thrusts[1,0], self.F_max)

        return motor_thrusts

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u