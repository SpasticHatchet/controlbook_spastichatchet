#%%

import sympy as sp

from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

#%%
# This is the full process through chapter 6 of design study B

t, m1, m2, ell, g = sp.symbols('t, m1, m2, ell, g')
theta, z, F = dynamicsymbols('theta, z, F')

q = sp.Matrix([[z],[theta]])
q_dot = q.diff(t)

half = sp.Rational(1, 2)

p1 = sp.Matrix([[z + half*ell*sp.sin(theta)], [half*ell*sp.cos(theta)], [0]])
v1 = p1.diff(t)

p2 = sp.Matrix([[z], [0], [0]])
v2 = p2.diff(t)

twelfth = sp.Rational(1, 12)

J = sp.diag(0, twelfth*(m1*ell**2), twelfth * (m1*ell**2))

omega_mat = sp.Matrix([[0],[0],[theta.diff(t)]])

K_1_trans = half * m1 * v1.T @ v1
K_1_rot = half * omega_mat.T @ J @ omega_mat
K_2 = half * m2 * v2.T @ v2

K = sp.simplify(K_1_trans + K_1_rot + K_2)
K = K[0,0]

display(Math(vlatex(K)))

x = sp.Matrix([[z], [theta], [z.diff(t)], [theta.diff(t)]])
x_dot = x.diff(t)


# %%
