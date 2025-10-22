#%%
import sympy as sp
from sympy import roots
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display


#%%
# from sympy import symbols, Eq, solve

# Define the variable and constants
s = sp.symbols('s')
b, m, k, k_d, k_p = sp.symbols('b m k k_d k_p')

params_dict = {
    b : 0.5,
    m : 5,
    k : 3
}

# Define the quadratic equation
b0 = 1/m
a1 = b/m
a0 = k/m

# D.7a
ol_char_eqn = sp.simplify(sp.Eq(s**2 + a1*s + a0, 0))

poles = sp.solve(ol_char_eqn, s)
poles_numerical = poles.subs[0](params_dict)
display(Math(vlatex(poles_numerical)))

#%%

cl_char_equation = sp.Eq(s**2 + (a1 + b0*k_d)*s + (a0 + b0*k_p), 0)


#solutions_params = cl_char_equation.subs(params_dict)

# Solve the equation for s
solutions = sp.solve(cl_char_equation, s)

#solutions = sp.simplify(solutions)

# Display the solutions
print('Closed-loop poles:')
display(Math(vlatex(solutions)))
# %%

# p1_eq = sp.Eq(solutions[0], sp.symbols('p_1'))
# p2_eq = sp.Eq(solutions[1], sp.symbols('p_2'))
p1_eq = sp.Eq(solutions[0], -1.5)
p2_eq = sp.Eq(solutions[1], -1)

#display(Math(vlatex(p1_eq)))
#display(Math(vlatex(p2_eq)))

#poles = sp.Matrix([[p1_eq],[p2_eq]])
#display(Math(vlatex(poles)))

# Solve the two scalar equations for k_p and k_d by passing them as a list
kp_kd = sp.solve([p1_eq, p2_eq], (k_p, k_d))

print('Kp, Kd:')
display(Math(vlatex(kp_kd)))

# Plugging in our values,
# k_p = 4.5, k_d = 12.0
# %%

# omega_nat = 2.2 / tr
# tr = 2 s
# ergo, omega_nat = 2.2/2

display(Math(vlatex(sp.symbols('omega_n'))))
print(2.2/2.0)
# %%
