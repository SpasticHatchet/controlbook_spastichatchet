#%%
import massParam as P
import numpy as np
import control as control

#%%
ctrlMat = control.ctrb(P.Amat, P.Bmat)

n = np.size(P.Amat, 1)

rank = np.linalg.matrix_rank(ctrlMat)

#%%
print(P.Bmat)
print(P.Amat @ P.Bmat)
print(ctrlMat)
print(n)
print(rank)
# %%

b = P.b
k = P.k
m = P.m

Kmat = np.matrix([[-0.122 -0.029108]])

A_BK = np.linalg.inv(P.Amat - (P.Bmat @ Kmat))

Kr = -1 / (P.Cmat @ A_BK @ P.Bmat)
print(Kr)
# %%
