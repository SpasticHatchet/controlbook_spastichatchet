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
t, m = symbols('t, m')
z = dynamicsymbols('z')

#defining generalized coords and derivatives
q = Matrix([z])
qdot = q.diff(t)

#defining the kinetic energy
p = Matrix([z])
v = p.diff(t)

half = sp.Rational(1,2)

K = simplify(half*m*v**2)

display(Math(vlatex(K)))