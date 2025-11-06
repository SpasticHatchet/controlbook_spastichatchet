#%%
import numpy as np
import VTOLParam as P
import control as control

#%%
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

C_AB_lon = control.ctrb(P.Amat_lon, P.Bmat_lon)

n_lon = np.size(P.Amat_lon, 1)
rank_lon = np.linalg.matrix_rank(C_AB_lon)
detC_AB_lon = np.linalg.det(C_AB_lon)
#%% 
# print(n_lon)
# print(rank_lon)
# print(detC_AB_lon)

K_lon = control.place(P.Amat_lon, P.Bmat_lon, des_poles_lon)

print("K_lon:")
print(K_lon)

A_BK = P.Amat_lon - P.Bmat_lon @ K_lon
# print("A-BK:")
# print(A_BK)

Krh = -1.0 / (P.Cmat_lon[0] @ np.linalg.inv(A_BK) @ P.Bmat_lon)

print("Krh:")
print(Krh[0,0])

# lateral
des_char_poly_lat = np.convolve([1, 2 * zeta_z * omegan_z, omegan_z**2],
                                [1, 2 * zeta_theta * omegan_theta, omegan_theta**2])
# print(des_char_poly_lat)
des_poles_lat = np.roots(des_char_poly_lat)
# print(des_poles_lat)

C_AB_lat = control.ctrb(P.Amat_lat, P.Bmat_lat)

n_lat = np.size(P.Amat_lat, 1)
rank_lat = np.linalg.matrix_rank(C_AB_lat)
detC_AB_lat = np.linalg.det(C_AB_lat)

# print(n_lat)
# print(rank_lat)
# print(detC_AB_lat)

K_lat = control.place(P.Amat_lat, P.Bmat_lat, des_poles_lat)

print("K_lat:")
print(K_lat)

A_BK_lat = P.Amat_lat - P.Bmat_lat @ K_lat

# print("A-BK_lat:")
# print(A_BK_lat)

Krz = -1.0 / (P.Cmat_lat[0] @ np.linalg.inv(A_BK_lat) @ P.Bmat_lat)

print("Krz:")
print(Krz[0,0])