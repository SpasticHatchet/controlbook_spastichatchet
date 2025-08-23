#%%
from hw02_arm_finding_KE import *

#%%[markdown]
# The code imported from above shows how we defined q, q_dot, and necessary system parameters. 
# Then we used position, velocity, and angular velocity to calculate kinetic energy.

#%%
#defining potential energy 
g = symbols('g')

P =  m*g*ell/2.0*sin(theta)   #this is "mgh", where "h" is a function of generalized coordinate "q"

# can also do the following to get the same answer
#   g_vec = Matrix([[0], [g], [0]])  # defining gravity in the direction that increases potential energy
#   p1 = Matrix([[ell/2*cos(theta)], [ell/2*sin(theta)], [0]])
#   P = m*g_vec.T@p1
#   P = P[0,0]

#calculate the lagrangian, using simplify intermittently can help the equations to be
#simpler, there are also options for factoring and grouping if you look at the sympy
#documentation.
L = simplify(K-P)

#%%
# Solution for Euler-Lagrange equations, but this does not include right-hand side (like friction and tau)
EL_case_studyA = simplify( diff(diff(L, qdot), t) - diff(L, q) )

display(Math(vlatex(EL_case_studyA)))



#%% 
############################################################
### Including friction and generalized forces, then solving for highest order derivatives
############################################################

# these are just convenience variables
thetad = theta.diff(t)
thetadd = thetad.diff(t)

# defining symbols for external force and friction
tau, b = symbols('tau, b')

# defining the right-hand side of the equation and combining it with E-L part
RHS = Matrix([[tau - b*thetad]])
full_eom = EL_case_studyA - RHS

# finding and assigning zdd and thetadd
# if our eom were more complicated, we could rearrange, solve for the mass matrix, and invert it to move it to the other side and find qdd and thetadd
result = simplify(sp.solve(full_eom, (thetadd)))

#TODO - add an example of finding the same thing, but not using sp.solve


# result is a Python dictionary, we get to the entries we are interested in
# by using the name of the variable that we were solving for
thetadd_eom = result[thetadd] # EOM for thetadd, as a function of states and inputs

display(Math(vlatex(thetadd_eom)))


#%% [markdown]
# OK, now we can get the state variable form of the equations of motion.

#%%
import armParam as P
import numpy as np

# defining fixed parameters that are not states or inputs (like g, ell, m, b)
# can be done like follows:
# params = [(m, P.m), (ell, P.ell), (g, P.g), (b, P.b)]  

# but in this example, I want to keep the masses, length, and damping as variables so
# that I can simulate uncertainty in those parameters in real life. 
params = [(g, P.g)]


# substituting parameters into the equations of motion
thetadd_eom = thetadd_eom.subs(params)

# now defining the state variables that will be passed into f(x,u) 
state = np.array([theta, thetad])
ctrl_input = np.array([tau])

# defining the function that will be called to get the derivatives of the states
state_dot = np.array([thetad, thetadd_eom])


#%%
import numpy as np

# converting the function to a callable function that uses numpy to evaluate and 
# return a list of state derivatives
eom = sp.lambdify([state, ctrl_input, m, ell, b], state_dot, 'numpy')

# calling the function as a test to see if it works:
cur_state = np.array([0, 0])
cur_input = np.array([1])
print("x_dot = ", eom(cur_state, cur_input, P.m, P.ell, P.b))




#%% [markdown] 
# The next step is to save this function "f" so that we can use it with a numerical integrator, like 
# scipy.integrate.ivp.solve_ivp or the rk4 functions in the case studies. To save this function, we do the
# following:

#%%
import dill   # we may have to install dill using pip "python -m pip install dill"
dill.settings['recurse'] = True
dill.dump(eom, open("eom_case_study_A", "wb"))   #takes the function name "eom" and writes it to a file called "eom_case_study_B"


# we can then reload the function and test it again to make sure it works: 
eom_test = dill.load(open("eom_case_study_A", "rb"))

print("x_dot after reloading = ", eom_test(cur_state, cur_input, P.m, P.ell, P.b))
