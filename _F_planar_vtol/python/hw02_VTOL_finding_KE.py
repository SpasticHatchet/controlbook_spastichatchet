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
t, ml, mr, mc, d, g, w = symbols('t, m_l, m_r, m_c, d, g, omega_c')
h = dynamicsymbols('h')
z = dynamicsymbols('z')
theta = dynamicsymbols('theta')

#defining generalized coords and derivatives
q = Matrix([[z], [theta], [0]])
qdot = q.diff(t)

#defining the kinetic energy
pc = Matrix([[z], [h], [0]])
pr = Matrix([[z + d*cos(theta)], [h + d*sin(theta)], [0]])
pl = Matrix([[z - d*cos(theta)], [h - d*sin(theta)], [0]])
vc = pc.diff(t)
vr = pr.diff(t)
vl = pl.diff(t)



half = sp.Rational(1,2)
#ml = mr
K = simplify(half*mr*vl.T@vl + half*mc*vc.T@vc + half*mr*vr.T@vr)

# just grabbing the scalar inside this matrix so that we can do L = K-P, since P is
K = K[0,0]

display(Math(vlatex(K)))