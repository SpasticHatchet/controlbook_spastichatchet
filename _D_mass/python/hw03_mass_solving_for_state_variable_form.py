from hw02_mass_finding_KE import *

k, F, b = symbols('k, F, b')

P = Matrix([[half*k*z**2]])

L = simplify(K - P)
L = L[0,0]

# display(Math(vlatex(L)))

EL_case_studyD = simplify( diff(diff(L, qdot), t) - diff(L, q) )
EL_case_studyD = EL_case_studyD[0,0]

# display(Math(vlatex(EL_case_studyD)))

tau = Matrix([[F]])

RHS = tau - b*qdot

# display(Math(vlatex(RHS)))

EOM = EL_case_studyD - RHS[0,0]

display(Math(vlatex(EOM)))