#%%
from h3_generate_E_L import * 
from IPython.display import display, Math
from sympy.physics.vector import vlatex

#%%
# TODO - calculate the full equations of motion for the hummingbird, you can do this
# as a single equation = 0 (i.e. LHS - RHS), or you can treat each part individually and
# do the calculations for LHS and RHS separately.

# We don't yet have a matrix for Beta, so let's create that first.
B = beta * sp.diag(1,1,1)

display(Math(vlatex(tau)))


RHS = tau - (B @ q_dot)
LHS = (M @ q_dot.diff(t)) + C + dP_dq

#display(Math(vlatex(RHS[1])))
#display(Math(vlatex(LHS[1])))

#%%[markdown]
# Longitudinal Dynamics are derived in this section: 

#%%
# step 0 - just print the equations for theta_ddot

display(Math(vlatex(sp.expand_trig(LHS[1]))))
display(Math(vlatex(sp.expand_trig(RHS[1]))))


#%%
# step 1 - substitute zero in for all the suggested variables
full_eom = sp.Eq(LHS[1], RHS[1])
#full_eom_oneside = LHS[1] - RHS[1]

zeroes_dict = {
    phi : 0,
    phi.diff(t) : 0,
    psi.diff(t) : 0,
    psi.diff(t, 2) : 0,
    J1x : 0 # J1x is being treated as time varying for some reason, so I'm just zeroing it out
}
zeroes_eom = full_eom.subs(zeroes_dict)
#zeroes_eom_os = full_eom_oneside.subs(zeroes_dict)

display('Longitudinal Dynamics (step 1):')
display(Math(vlatex(zeroes_eom)))

#%%
# step 2 - substitute F and Tau for f_l, f_r, and d*(f_l-f_r)
F, tau2 = sp.symbols('F, tau')

zeroes_eom_F = zeroes_eom.subs(f_l + f_r, F)
zeroes_eom_d = zeroes_eom_F.subs(d*(f_l - f_r), tau2)
#zeroes_eom_F_os = zeroes_eom_os.subs(f_l + f_r, F)
#zeroes_eom_d_os = zeroes_eom_F_os.subs(d*(f_l - f_r), tau2)

display('Longitudinal Dynamics (step 2):')
display(Math(vlatex(zeroes_eom_d)))

#%%
# step 3 - find the feedback linearization force F_fl
display('Feedback Linearization Force (step 3):')

linear_zeroes = {
    theta.diff(t) : 0,
    theta.diff(t,2) : 0,
}

nonlinear = zeroes_eom_d.subs(linear_zeroes)

F_fl = sp.solve(nonlinear, F)[0]
display(Math(vlatex(sp.Eq(sp.symbols('F_fl'), F_fl))))
#display(Math(vlatex(nonlinear)))

#%%
# step 4 - find the linearized equation of motion for theta_ddot
display('Linearized Longitudinal Dynamics (step 4):')

Fc = sp.symbols('F_control')
linearized = zeroes_eom_d.subs(F, Fc + F_fl)
oneside = sp.simplify(linearized.lhs - linearized.rhs)
result = sp.collect(oneside, theta.diff(t,2))

solved = sp.solve(result, theta.diff(t,2))[0]
#display(Math(vlatex(sp.simplify(solved))))

# Set friction coefficient b to 0, which makes the B matrix 0
theta_ddot_rhs = solved.subs(beta, 0)
theta_ddot_eq = sp.Eq(theta.diff(t,2), theta_ddot_rhs)
display(Math(vlatex(theta_ddot_eq)))

#%% [markdown]
# Lateral Dynamics are derived in this section:

#%%
# step 1 - substitute zero in for all the suggested variables
display('Lateral Dynamics (step 1):')
#display(Math(vlatex()))
#display(Math(vlatex(sp.simplify(LHS))))
#display(Math(vlatex(sp.simplify(RHS))))

full_eom_r0 = sp.Eq(LHS[0], RHS[0])
full_eom_r2 = sp.Eq(LHS[2], RHS[2])


zeroes_dict_lat = {
    theta : 0,
    theta.diff(t) : 0,
    theta.diff(t, 2) : 0
}

zeroes_eom_r0 = full_eom_r0.subs(zeroes_dict_lat)
zeroes_eom_r2 = full_eom_r2.subs(zeroes_dict_lat)

#Make beta 0, replace fl+fr with F, replace d(fl-fr) with tau2
subs_dict_lat = {
    f_l + f_r : F,
    d*(f_l - f_r) : tau2,
    beta: 0
}
replacement_r0 = sp.trigsimp(zeroes_eom_r0.subs(subs_dict_lat))
replacement_r2 = sp.trigsimp(zeroes_eom_r2.subs(subs_dict_lat))

display(Math(vlatex(replacement_r0)))
display(Math(vlatex(replacement_r2)))

#%%
# step 2 - find the equilibrium force F_e
zero_r1_thetas = zeroes_eom_d.subs(zeroes_dict_lat)
F_e = sp.solve(zero_r1_thetas, F)[0]

display(Math(vlatex(sp.Eq(sp.symbols('F_e'), F_e))))


#%%
# step 3 - find the equilibrium variables Tau_e and phi_e

#display(Math(vlatex()))
# by inspection Tau_e = 0, phi_e = 0 or pi, but 0 makes the most sense for the hb

# To solve for Tau_e and phi_e, we set phi_ddot, psi_ddot, phi_dot, and psi_dot to zero
equilibrium_zeroes = {
    phi.diff(t) : 0,
    phi.diff(t, 2) : 0,
    psi.diff(t) : 0,
    psi.diff(t, 2) : 0

}
# Sub those into our lateral equations of motion
zero_r0_eq = replacement_r0.subs(equilibrium_zeroes)
zero_r2_eq = replacement_r2.subs(equilibrium_zeroes)

# Solve for Tau_e and phi_e by plugging in F_e
zero_r0_eq = zero_r0_eq.subs(F, F_e)
zero_r2_eq = zero_r2_eq.subs(F, F_e)
Tau_e = sp.solve(zero_r0_eq, tau2)[0]
phi_e = sp.solve(zero_r2_eq, phi)[0]

display(Math(vlatex(sp.Eq(sp.symbols('Tau_e'), Tau_e))))
display(Math(vlatex(sp.Eq(sp.symbols('phi_e'), phi_e))))

#%%
# step 4 - find the linearized equation of motion for phi_ddot and psi_ddot

# Define x and u vectors
x_lat = sp.Matrix([[psi], [phi], [psi.diff(t)], [phi.diff(t)]])
u_lat = sp.Matrix([tau2])

phi_ddot_eom = sp.solve(replacement_r0, phi.diff(t,2))[0]
psi_ddot_eom = sp.solve(replacement_r2, psi.diff(t,2))[0]

x_dot = sp.Matrix([psi.diff(t), phi.diff(t), psi_ddot_eom, phi_ddot_eom])

A_mat_pre = x_dot.jacobian(x_lat)
B_mat_pre = x_dot.jacobian(u_lat)

x_lat_eq = {
    psi : 0,
    phi : phi_e,
    psi.diff(t) : 0,
    phi.diff(t) : 0,
    F : sp.symbols('F_e')
}
u_lat_eq = {
    tau2 : Tau_e
}

A_mat = A_mat_pre.subs(x_lat_eq | u_lat_eq)
B_mat = B_mat_pre.subs(x_lat_eq | u_lat_eq)

display('Linearized Lateral Dynamics (step 4): ')
display(Math(vlatex(sp.simplify(A_mat))))
display(Math(vlatex(sp.simplify(B_mat))))

# HINT: I would suggest finding equations for phi_ddot and psi_ddot, then 
# you can linearize them like we have shown in the case studies using the 
# "jacobian" function, and substituting equilibrium values in. 



# %%
