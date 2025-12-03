import numpy as np
import VTOLParam as P
from scipy import signal
import control as cnt

class ctrlLQR:
    def __init__(self):
        # tuning parameters
        Q_lon = np.diag([0.01, 100.0, 0.01, 5000.0, 1000.0])
        R_lon = np.diag([0.5])  # [tau]
        Q_lat = np.diag([30.0, 100.0, 1.0, 1.0, 1000.0])
        R_lat = np.diag([1.0])  # [F]
        # tr_h    = 5.0 
        # zeta_h  = 0.95
        # wn_h    = 0.5*np.pi/(tr_h*np.sqrt(1-zeta_h**2))

        # tr_z    = 5.0 
        # zeta_z  = 0.95
        # wn_z    = 0.5*np.pi/(tr_z*np.sqrt(1-zeta_z**2))

        # tr_th   = tr_z/10.0 
        # zeta_th = 0.95
        # wn_th   = 0.5*np.pi/(tr_th*np.sqrt(1-zeta_th**2))

        # integrator_h = 1.0
        # integrator_z = 1.0

        # # observer gains
        # wn_h_obs    = 10.0*wn_h
        # wn_z_obs    = 10.0*wn_z
        # wn_th_obs   = 10.0*wn_th

        # State Space Equations
        self.Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force 
        self.A_lon = np.array([[0.0, 1.0],
                        [0.0, 0.0]])
        self.B_lon = np.array([[0.0],
                        [1.0 / (P.mc + 2.0 * P.mr)]])
        self.C_lon = np.array([[1.0, 0.0]])

        self.A_lat = np.array([[0.0, 0.0, 1.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0],
                        [0.0, -self.Fe / (P.mc + 2.0 * P.mr), -(P.mu / (P.mc + 2.0 * P.mr)), 0.0],
                        [0.0, 0.0, 0.0, 0.0]])
        self.B_lat = np.array([[0.0],
                        [0.0],
                        [0.0],
                        [1.0 / (P.Jc + 2 * P.mr * P.d ** 2)]])
        self.C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0]])
        
        # form augmented system
        #   for longitudinal system
        Cr_lon = self.C_lon
        A1_lon = np.zeros((self.A_lon.shape[0]+1, self.A_lon.shape[0]+1))
        A1_lon[0:2, 0:2] = self.A_lon
        A1_lon[2, 0:2] = -Cr_lon

        B1_lon = np.zeros((self.B_lon.shape[0]+1, 1))        
        B1_lon[0:-1, :] = self.B_lon

        #   for lateral system
        Cr_lat = self.C_lat[0,:]
        A1_lat = np.zeros((self.A_lat.shape[0]+1, self.A_lat.shape[0]+1))
        A1_lat[0:4, 0:4] = self.A_lat
        A1_lat[4, 0:4] = -Cr_lat
        B1_lat = np.zeros((self.B_lat.shape[0]+1, 1))
        B1_lat[0:-1, :] = self.B_lat

        # gain calculation
        # des_char_poly_lon = np.convolve([1.0, 2.0 * zeta_h * wn_h, wn_h ** 2],
        #                                 [1, integrator_h])
        # des_poles_lon = np.roots(des_char_poly_lon)
        
        # des_char_poly_lat = np.convolve(
        #     np.convolve([1.0, 2.0 * zeta_z * wn_z, wn_z ** 2],
        #                 [1.0, 2.0 * zeta_th * wn_th, wn_th ** 2]),
        #     [1, integrator_z])
        # des_poles_lat = np.roots(des_char_poly_lat)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != A1_lon.shape[0]:
            print("The longitudinal system is not controllable")
        else:
            K1_lon_lqr = cnt.lqr(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]

        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != A1_lat.shape[0]:
            print("The lateral system is not controllable")
        else:
            K1_lat = cnt.place(A1_lat, B1_lat, des_poles_lat)
            self.K_lat = K1_lat[0][0:4]
            self.ki_lat = K1_lat[0][4]

        # compute observer poles
        des_obs_poles_lon = np.roots([1.0, 2.0*zeta_h*wn_h_obs, wn_h_obs**2])
        des_obs_poles_lat = np.convolve([1.0, 2.0*zeta_z*wn_z_obs, wn_z_obs**2],
                                        [1.0, 2.0*zeta_th*wn_th_obs, wn_th_obs**2])
        des_obs_poles_lat = np.roots(des_obs_poles_lat)

        # check if system is observable and calculate L for lon. and lat. systems
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lon.T, self.C_lon.T)) != self.A_lon.shape[0]:
            print("The longitudinal system is not observable")
        else:
            self.L_lon = signal.place_poles(self.A_lon.T, self.C_lon.T, des_obs_poles_lon).gain_matrix.T

        if np.linalg.matrix_rank(cnt.ctrb(self.A_lat.T, self.C_lat.T)) != self.A_lat.shape[0]:
            print("The lateral system is not observable")
        else:
            self.L_lat = signal.place_poles(self.A_lat.T, self.C_lat.T, des_obs_poles_lat).gain_matrix.T

        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)
        print('L_lon^T: ', self.L_lon.T)
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        print('L_lat^T: ', self.L_lat.T)

        # initialize estimated states
        self.xhat_lon = np.array([[0.0], [0.0]])
        self.xhat_lat = np.array([[0.0], [0.0], [0.0], [0.0]])

        self.integrator_z = 0.0  # integrator on position z
        self.error_z_prev = 0.0  # error signal delayed by 1 sample
        self.integrator_h = 0.0  # integrator on altitude h
        self.error_h_prev = 0.0  # error signal delayed by 1 sample
        self.F_prev = 0.0  # Force signal delayed by 1 sample
        self.tau_prev = 0.0  # torque signal delayed by 1 sample
        self.F_limit = P.F_max * 2.0
        self.tau_limit = P.F_max * P.d * 2.0
        self.Ts = P.Ts

    def update(self, r, y):
        z_r = r[0][0]
        h_r = r[1][0]

        # assigning measurements
        y_lat = np.array([[y[0][0]],[y[2][0]]])
        y_lon = y[1][0]

        # update the observers to control with the estimated states 
        xhat_lat = self.update_lat_observer(y_lat)
        xhat_lon = self.update_lon_observer(y_lon)
        z_hat = xhat_lat[0][0]
        h_hat = xhat_lon[0][0]
        theta_hat = xhat_lat[1][0]

        # integrate error
        error_z = z_r - z_hat
        self.integrator_z += (P.Ts/2.0)*(error_z + self.error_z_prev)
        self.error_z_prev = error_z

        error_h = h_r - h_hat
        self.integrator_h += (P.Ts/2.0)*(error_h + self.error_h_prev)
        self.error_h_prev = error_h

        # Compute the state feedback controllers
        F_tilde = -self.K_lon @ xhat_lon - self.ki_lon * self.integrator_h
        F = self.Fe + F_tilde[0]
        F_sat = saturate(F, self.F_limit)
        self.integratorAntiWindup(F_sat, F, self.ki_lon, self.integrator_h)

        tau = -self.K_lat @ xhat_lat - self.ki_lat*self.integrator_z
        tau_sat = saturate(tau[0], self.tau_limit)
        self.integratorAntiWindup(tau_sat, tau, self.ki_lat, self.integrator_z)

        u = np.array([[F_sat], [tau_sat]])
        self.F_prev = F_sat
        self.tau_prev = tau_sat

        return u, xhat_lat, xhat_lon
    
    def integratorAntiWindup(self, u_sat, u_unsat, ki, integrator):
        if ki != 0.0:
            integrator = integrator + P.Ts/ki*(u_sat-u_unsat)

    def update_lat_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lat(self.xhat_lat, y_m)
        F2 = self.observer_f_lat(self.xhat_lat + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lat(self.xhat_lat + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lat(self.xhat_lat + self.Ts * F3, y_m)
        self.xhat_lat += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        return self.xhat_lat

    def observer_f_lat(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xhat_dot = self.A_lat @ x_hat \
                   + self.B_lat * self.tau_prev \
                   + self.L_lat @ (y_m - self.C_lat @ x_hat)
        
        return xhat_dot

    def update_lon_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lon(self.xhat_lon, y_m)
        F2 = self.observer_f_lon(self.xhat_lon + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lon(self.xhat_lon + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lon(self.xhat_lon + self.Ts * F3, y_m)
        self.xhat_lon += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        return self.xhat_lon

    def observer_f_lon(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xhat_dot = self.A_lon @ x_hat \
                   + self.B_lon * (self.F_prev - self.Fe) \
                   + self.L_lon @ (y_m - self.C_lon @ x_hat)
        
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit*np.sign(u)
    return u

