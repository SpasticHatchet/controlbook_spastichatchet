import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlPD import ctrlPD

# instantiate blockbeam, controller, and reference classes
blockbeam = blockbeamDynamics(alpha=0.0)
ctrl = ctrlPD()
reference = signalGenerator(amplitude=0.15, frequency=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
y = blockbeam.h()  # output of system at start of simulation

# for part e), we can uncomment below
# blockbeam.state[1,0] = 10.0*np.pi/180.0

while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        
        r = reference.square(t)  + 0.25 # reference input
        x = blockbeam.state  # use state instead of output
        u = ctrl.update(r, x)  # update controller
        y = blockbeam.update(u)  # propagate system
        t += P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, blockbeam.state, u, r)
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
