import numpy as np
import massParam as P

class massController_PD:
    def __init__(self):
        self.Kp = -P.k + 1.5*P.m
        self.Kd = -P.b + 2.5*P.m
        
    def update(self, r, state):

        z = state[0,0]
        zdot = state[1,0]

        u = (self.Kp * (r - z)) - (self.Kd * zdot)

        return u
