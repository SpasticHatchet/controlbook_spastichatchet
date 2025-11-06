import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedback:
    def __init__(self):

        self.F_max = P.F_max

        tr_theta = 0.8
        omegan_theta = 2.2 / tr_theta
        zeta_theta = 0.707

        tr_z = 10 * tr_theta
        omegan_z = 2.2 / tr_z
        zeta_z = 0.707

        tr_h = 8.0
        omegan_h = 2.2 / tr_h
        zeta_h = 0.707

        # longitudinal
        des_char_poly_lon = [1, 2 * zeta_h * omegan_h, omegan_h**2]

        # print(des_char_poly_lon)

        des_poles_lon = np.roots(des_char_poly_lon)

        # print(des_poles_lon)

        C_AB_lon = cnt.ctrb(P.Amat_lon, P.Bmat_lon)

        n_lon = np.size(P.Amat_lon, 1)
        rank_lon = np.linalg.matrix_rank(C_AB_lon)
        detC_AB_lon = np.linalg.det(C_AB_lon)
        #%% 
        # print(n_lon)
        # print(rank_lon)
        # print(detC_AB_lon)

        self.K_lon = cnt.place(P.Amat_lon, P.Bmat_lon, des_poles_lon)

        print("K_lon:")
        print(self.K_lon)

        A_BK = P.Amat_lon - P.Bmat_lon @ self.K_lon
        # print("A-BK:")
        # print(A_BK)

        self.Krh = -1.0 / (P.Cmat_lon[0] @ np.linalg.inv(A_BK) @ P.Bmat_lon)

        print("Krh:")
        print(self.Krh[0,0])

        # lateral
        des_char_poly_lat = np.convolve([1, 2 * zeta_z * omegan_z, omegan_z**2],
                                        [1, 2 * zeta_theta * omegan_theta, omegan_theta**2])
        # print(des_char_poly_lat)
        des_poles_lat = np.roots(des_char_poly_lat)
        # print(des_poles_lat)

        C_AB_lat = cnt.ctrb(P.Amat_lat, P.Bmat_lat)

        n_lat = np.size(P.Amat_lat, 1)
        rank_lat = np.linalg.matrix_rank(C_AB_lat)
        detC_AB_lat = np.linalg.det(C_AB_lat)

        # print(n_lat)
        # print(rank_lat)
        # print(detC_AB_lat)

        self.K_lat = cnt.place(P.Amat_lat, P.Bmat_lat, des_poles_lat)

        print("K_lat:")
        print(self.K_lat)

        A_BK_lat = P.Amat_lat - P.Bmat_lat @ self.K_lat

        # print("A-BK_lat:")
        # print(A_BK_lat)

        self.Krz = -1.0 / (P.Cmat_lat[0] @ np.linalg.inv(A_BK_lat) @ P.Bmat_lat)

        print("Krz:")
        print(self.Krz[0,0])

    def update(self, r, x):
        
        z_r = r[0]
        h_r = r[1]
        
        # Extract individual states for readability (using .item() for scalars from NumPy matrix)
        z = x[0]     # x: z_tilde
        h = x[1]     # x[1]: h_tilde
        theta = x[2] # x[2]: theta_tilde
        zdot = x[3]  # x[3]: z_dot_tilde
        hdot = x[4]  # x[4]: h_dot_tilde
        thetadot = x[5] # x[5]: theta_dot_tilde

        x_lon = np.array([h, hdot])

        print("x_lon:")
        print(x_lon)
        
        x_lat = np.array([[z], [theta], [zdot], [thetadot]]) 
        
        F_tilde_matrix = -self.K_lon @ x_lon + self.Krh * h_r
        
        tau_tilde_matrix = -self.K_lat @ x_lat + self.Krz * z_r

        F_total = P.Fe + F_tilde_matrix
        
        d_term = 1.0 / (2.0 * P.d) * tau_tilde_matrix
        
        fr = 0.5 * F_total + d_term
        fl = 0.5 * F_total - d_term

        fr_sat = self.saturate(fr, self.F_max)
        fl_sat = self.saturate(fl, self.F_max)

        return np.array([[fr_sat], [fl_sat]])

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u

