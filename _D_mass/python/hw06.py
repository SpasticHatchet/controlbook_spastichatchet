#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display


#%%
# from sympy import symbols, Eq, solve

# Define the variable and constants
s = sp.symbols('s')
b, m, k, k_d, k_p = sp.symbols('b m k k_d k_p')

# Define the quadratic equation
cl_char_equation = sp.Eq(s**2 + ((b/m) + (1/m)*k_d)*s + ((k/m) + (1/m)*k_p), 0)

# Solve the equation for s
solutions = sp.solve(cl_char_equation, s)
#solutions = sp.simplify(solutions)

# Display the solutions
print('Closed-loop poles:')
display(Math(vlatex(solutions)))
# %%

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
