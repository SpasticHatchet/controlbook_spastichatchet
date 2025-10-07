# %%
# these are libraries and functions that will help us with the calculations for h2
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display


def printsym(expr):
    """
    Renders latex from sympy expression when in a Jupyter-like environment. Does
    not print anything useful when running in a terminal. When in a terminal, use
    sp.pretty_print() instead.
    """
    display(Math(vlatex(expr)))


# %%
# These functions are defined in our helper functions file in the public repository
from EL_helper_functions import rotx, roty, rotz, calc_omega, find_coeffs

# Defining necessary symbols and variables to use in the calculations
t, ell_1, ell_2, ell_3x, ell_3y, ell_3z, m1, m2, m3 = sp.symbols(
    "t, ell_1, ell_2, ell_3x, ell_3y, ell_3z, m1, m2, m3"
)

# We have to tell SymPy that these are functions of time
theta = dynamicsymbols("theta")
psi = dynamicsymbols("psi")
phi = dynamicsymbols("phi")

q = sp.Matrix([[phi], [theta], [psi]])
q_dot = q.diff(t)


# %%
# TODO: define the position of each mass in body frame, then rotate
# it into the world or inertial frame.

p1_in_b = sp.Matrix([[ell_1],[0],[0]])
p1_in_w = sp.Matrix([[ell_1*sp.cos(theta)*sp.cos(psi)],[ell_1*sp.cos(theta)*sp.sin(psi)],[-ell_1*sp.sin(theta)]])

p2_in_2 = sp.Matrix([[ell_2],[0],[0]])
p2_in_w = rotz(psi)@roty(theta)@p2_in_2

p3_in_1 = sp.Matrix([[ell_3x], [ell_3y], [ell_3z]])
p3_in_w = rotz(psi)@p3_in_1


# %%
# TODO: take the time derivative of the position vectors to get the linear velocity,
# then use the "find_coeffs" function to calculate the "V_i" matrices

v1_in_w = sp.diff(p1_in_w)
V1 = find_coeffs(v1_in_w, q_dot)

v2_in_w = sp.diff(p2_in_w)
V2 = find_coeffs(v2_in_w, q_dot)

v3_in_w = sp.diff(p3_in_w)
V3 = find_coeffs(v3_in_w, q_dot)


# %%
# TODO: use the "rotx", "roty", and "rotz" functions to calculate the rotation matrices
# for each rigid body

R1 =  rotz(psi)@roty(theta)@rotx(phi)
R2 =  rotz(psi)@roty(theta)
R3 =  rotz(psi)


# %%
# TODO: use the "calc_omega" function with the rotation matrices to calculate the
# angular velocity of each rigid body

omega_1 = calc_omega(R1)
omega_2 = calc_omega(R2)
omega_3 = calc_omega(R3)


# %%
# TODO: use the "find_coeffs" function to calculate the "W_i" matrices

W1 = find_coeffs(omega_1, q_dot)
W2 = find_coeffs(omega_2, q_dot)
W3 = find_coeffs(omega_3, q_dot)


# %%
# TODO: define the inertia tensors for each rigid body

J1z, J1y, J1x, J2z, J2y, J2x, J3z, J3y, J3x = sp.symbols('J1_z, J1_y, J1_x, J2_z, J2_y, J2_x, J3_z, J3_y, J3_x')

J1 = sp.diag(J1x, J1y, J1z)

J2 = sp.diag(J2x, J2y, J2z)

J3 = sp.diag(J3x, J3y, J3z)


# %%
# TODO: calculate M using the masses and the V, W, R, and J matrices

# M = sp.zeros(3, 3)
M = (m1*V1.T@V1 + W1.T@R1@J1@R1.T@W1 + 
     m2*V2.T@V2 + W2.T@R2@J2@R2.T@W2 + 
     m3*V3.T@V3 + W3.T@R3@J3@R3.T@W3)
# M = M +


# %%
# Simplify and display the result
M = sp.trigsimp(M)
printsym(M)

# %%

