#%%
import numpy as np
import blockbeamParam as P
import control as control

#%%
tr_theta = 0.25
tr_z = 10 * tr_theta
omegan_z = 2.2 / tr_z # 2.2/tr
zeta_z = 0.707
omegan_theta = 2.2 / tr_theta
zeta_theta = 0.707

des_char_poly = np.convolve(
    [1, 2 * zeta_z * omegan_z, omegan_z**2],
    [1, 2 * zeta_theta * omegan_theta, omegan_theta**2])


print(des_char_poly)


des_poles = np.roots(des_char_poly)


print(des_poles)

C_AB = control.ctrb(P.Amat, P.Bmat)

n = np.size(P.Amat, 1)
rank = np.linalg.matrix_rank(C_AB)
detC_AB = np.linalg.det(C_AB)
#%% 
# print(n)
# print(rank)
# print(detC_AB)

K = control.place(P.Amat, P.Bmat, des_poles)


print(K)

A_BK = P.Amat - P.Bmat @ K
print("A-BK:")
print(A_BK)

Kr = -1.0 / (P.Cmat[0] @ np.linalg.inv(A_BK) @ P.Bmat)

print("Kr:")
print(Kr)
# %%
