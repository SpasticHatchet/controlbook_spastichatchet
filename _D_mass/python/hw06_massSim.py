import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from massController_PD import massController_PD

massCtrl = massController_PD()
massDyn = massDynamics()
z_disp = signalGenerator(amplitude = 1, frequency = 0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start # time starts at t_start
while t < P.t_end: # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot: # updates control and dynamics at faster simulation rate
        
        r = z_disp.square(t)
        u = massCtrl.update(r, massDyn.state)
        y = massDyn.update(u) # Propagate the dynamics

        t += P.Ts # advance time by Ts
    
    # update animation and data plots
    animation.update(massDyn.state)
    dataPlot.update(t, massDyn.state, u, r)
    plt.pause(0.01) # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()