# import matplotlib.pyplot as plt
# import numpy as np
# import massParam as P
# from signalGenerator import signalGenerator
# from massAnimation import massAnimation
# from dataPlotter import dataPlotter

# # instantiate reference input classes, these are not actual values,
# # just values to allow us to plot
# z_plot = signalGenerator(amplitude=2.0*np.pi, frequency=0.1)
# # tau_plot = signalGenerator(amplitude=5, frequency=.5)

# # instantiate the simulation plots and animation
# dataPlot = dataPlotter()
# animation = massAnimation()

# t = P.t_start  # time starts at t_start
# while t < P.t_end:  # main simulation loop
#     # set variables
#     z = z_plot.sin(t)
#     # tau = tau_plot.sawtooth(t)

#     # update animation
#     z_dot = 0.0

#     # state is made of theta, and theta_dot
#     state = np.array([[z], [z_dot]])  

#     # update animation and plot 
#     animation.update(state)
#     dataPlot.update(t, state, 1.0)# tau)

#     # advance time by t_plot
#     t += P.t_plot
#     plt.pause(0.02)  # allow time for animation to draw

# # Keeps the program from closing until the user presses a button.
# print('Press key to close')
# plt.waitforbuttonpress()
# plt.close()

import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

mass = massDynamics()
force = signalGenerator(amplitude = 10.0, frequency = 1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start # time starts at t_start
while t < P.t_end: # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot: # updates control and dynamics at faster simulation rate
        u = force.sin(t)
        y = mass.update(u) # Propagate the dynamics
        t += P.Ts # advance time by Ts
    
    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, mass.state, u)
    plt.pause(0.01) # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()