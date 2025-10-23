#%%
from hw6_blockbeam_linearization import *
from sympy import eye, zeros, simplify

# we can also transform this to a transfer function if we define C and D matrices
C = Matrix([[1.0, 0, 0, 0], [0, 1.0, 0, 0]])
D = Matrix([[0], [0]])

# these are the two transfer functions for z and theta with respect to input F
s = symbols('s')

# this is the C(sI-A)^(-1)B + D  equation shown in class to find the transfer function
transfer_func = simplify(C@(s*eye(4)-A_lin).inv()@B_lin+D)   
print("Tranfer functions (without simplifying assumption):")
display(Math(vlatex(transfer_func)))

# now setting the m1*g term equal to zero as described in HW E.5
A_E5 = A_lin.subs([(m1*g, 0)])
B_E5 = B_lin.subs([(m1*g, 0)])

transfer_func_partC = simplify(C*(s*eye(4)-A_E5).inv()*B_E5+D)
print("Tranfer functions (with simplifying assumption):")
display(Math(vlatex(transfer_func_partC)))


# %%
