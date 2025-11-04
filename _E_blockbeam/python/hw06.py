#%%
import blockbeamParam as P

tr_th = 1.0       # Rise time for inner loop (theta)
zeta_th = 0.707       # inner loop Damping Coefficient

A = ((P.m2 * P.ell**2.0) / 3.0) + (P.m1 * (P.ell / 2.0)**2.0)
#---------------------------------------------------
#                    Inner Loop
#---------------------------------------------------
# parameters of the open loop transfer function
b0_th = P.ell / A
a0_th = 0.0
a1_th = 0.0

# coefficients for desired inner loop
wn_th = 2.2 / tr_th     # Natural frequency

# compute gains
kp_th = (wn_th**2 + a0_th) / b0_th
kd_th = (2.0 * zeta_th * wn_th) / b0_th

print(kp_th)
print(kd_th)
# %%
DCgain = b0_th * kp_th / (b0_th * kp_th + a0_th)