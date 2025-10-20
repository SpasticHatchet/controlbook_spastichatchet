# #%%

# the following commented code works correctly, but the book solution made me think I had done 
# something wrong. It was really just a difference in how things were simplified.

# import sympy as sp
# from sympy.physics.vector import dynamicsymbols
# from sympy.physics.vector.printing import vlatex
# from IPython.display import Math, display

# # importing these functions directly to make life a little easier and code a little more readable
# from sympy import sin, cos, Matrix, symbols, simplify, diff
# # init_printing(use_latex=True)


# #defining mathematical variables (called symbols in sp) and time varying functions like z and theta
# t, m1, m2, L, g = symbols('t, m1, m2, ell, g')

# # force = signalGenerator(amplitude=0.5, frequency=1.0, y_offset=11.5)
# # u = force.sin(t)

# #defining generalized coord and their derivatives
# theta = dynamicsymbols('theta')
# z = dynamicsymbols('z')
# q = Matrix([[z, theta]])
# qdot = q.diff(t)

# #to find KE, start with the position of the mass in the inertial frame, then find the velocity
# half = sp.Rational(1, 2)
# p1 = Matrix([[z*cos(theta)], [z*sin(theta)], [0]])
# p2 = Matrix([[half*L*cos(theta)], [half*L*sin(theta)], [0]])
# v1 = p1.diff(t)
# v2 = p2.diff(t)

# # define the rotation matrix describing the orientation of the rod relative to the inertial frame
# # R = Matrix([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])
# omega = Matrix([[0], [0], [sp.diff(theta)]])


# # find the angular velocity of the rod
# # Omeg_mat = simplify(R.diff(t)*R.T)

# # omega = Matrix([[Omeg_mat[2,1]], [Omeg_mat[0,2]], [Omeg_mat[1,0]]])
# # the following function does the same thing as above 
# # (extracting omega - vector - from Omeg_mat - a skew symmetric matrix)
# #   omega = Omeg_mat.vee()

# # next we define the inertia tensor for the rod
# twelfth = sp.Rational(1, 12)
# J = sp.diag(0, 0, twelfth*m2*L**2)

# # calculate the kinetic energy and display it
# K = simplify((half*m1*v1.T @ v1) + (half*m2*v2.T @ v2) + (half*omega.T @ J @ omega))
# #K = simplify(0.5*m*v.T @ v + 0.5*omega.T @ R @ J @ R.T @ omega)
# K = K[0,0]

# display(Math(vlatex(K)))

import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

###################################################################################
# find Kinetic Energy (with vectors/matrices) for case study E blockbeam
###################################################################################
# importing these functions directly to make life a little easier and code a little
from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

#defining mathematical variables (called symbols in sp) and time varying functions
t, m1, m2, ell, g = symbols('t, m1, m2, ell, g')
z = dynamicsymbols('z')
theta = dynamicsymbols('theta')

#defining generalized coords and derivatives
q = Matrix([[z], [theta]])
qdot = q.diff(t)

#defining the kinetic energy
p1 = Matrix([[z*cos(theta)], [z*sin(theta)], [0]])
p2 = Matrix([[ell/2*cos(theta)], [ell/2*sin(theta)], [0]])
v1 = p1.diff(t)
v2 = p2.diff(t)

# by inspection, we can find the angular velocity of the rod
omega = Matrix([[0], [0], [theta.diff(t)]])

# if we are uncomfortable with the above, we can use the rotation matrix describing
# the orientation of the beam relative to the inertial frame to find the angular ve
R = Matrix([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])
# Omeg_mat = R.diff(t)*R.T
# omega = Omeg_mat.vee()

# next we define the inertia tensor for the beam, modeled as a thin rod
#J = diag(0, m2*ell**2/12.0, m2*ell**2/12.0)
twelfth = sp.Rational(1, 12)
J = sp.diag(0, twelfth*m2*ell**2, twelfth*m2*ell**2)

half = sp.Rational(1,2)
K = simplify(half*m1*v1.T@v1 + half*m2*v2.T@v2 + half*omega.T@R@J@R.T@omega)

# just grabbing the scalar inside this matrix so that we can do L = K-P, since P is
K = K[0,0]

display(Math(vlatex(K)))