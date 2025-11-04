#%%
import sympy as sp
from sympy import roots
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

#%%

s = sp.symbols('s')
b, m, k, k_d, k_p = sp.symbols('b m k k_d k_p')

params_dict = {
    b : 0.5,
    m : 5,
    k : 3
}

b0 = 1/m
a1 = b/m
a0 = k/m

#%%

#D.7a
print('Open-loop poles:')
ol_char_eqn = sp.simplify(sp.Eq(s**2 + a1*s + a0, 0))
ol_poles = sp.solve(ol_char_eqn, s)

ol_p1 = ol_poles[0].subs(params_dict)
ol_p2 = ol_poles[1].subs(params_dict)

display(Math(vlatex(ol_poles)))
display(Math(vlatex([ol_p1, ol_p2])))

#%%

#D.7b
cl_char_eqn = sp.Eq(s**2 + (a1 + b0*k_d)*s + (a0 + b0*k_p), 0)
cl_poles = sp.solve(cl_char_eqn, s)

print('Closed-loop poles:')
cl_p1 = sp.simplify(cl_poles[0].subs(params_dict))
cl_p2 = sp.simplify(cl_poles[1].subs(params_dict))

display(Math(vlatex(cl_poles)))
display(Math(vlatex([cl_p1, cl_p2])))
# %%

#D.7c
p1_eq = sp.Eq(cl_poles[0], -1.5)
p2_eq = sp.Eq(cl_poles[1], -1)

kp_kd = sp.solve([p1_eq, p2_eq], (k_p, k_d))

# kp_num = kp_kd[0].subs(params_dict)
# kd_num = kp_kd[1].subs(params_dict)
# Plugging in our values, k_p = 4.5, k_d = 12.0

kp_num = 4.5
kd_num = 12.0

print('Kp, Kd:')
display(Math(vlatex(kp_kd)))
print(kp_num, kd_num)

# Run hw06_massSim.py for D.7d
