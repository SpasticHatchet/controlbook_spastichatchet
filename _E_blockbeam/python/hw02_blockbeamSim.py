# import matplotlib.pyplot as plt
# import numpy as np
# import blockbeamParam as P
# from signalGenerator import signalGenerator
# from blockbeamAnimation import blockbeamAnimation
# from dataPlotter import dataPlotter

# # instantiate reference input classes
# z_plot = signalGenerator(amplitude=0.2, frequency=0.1, y_offset=0.25)
# theta_plot = signalGenerator(amplitude=.05*np.pi, frequency=.5)
# # f_fakeValueGenerator = signalGenerator(amplitude=5, frequency=.5)

# # instantiate the simulation plots and animation
# dataPlot = dataPlotter()
# animation = blockbeamAnimation()

# t = P.t_start  # time starts at t_start
# while t < P.t_end:  # main simulation loop
#     # set variables
#     z = z_plot.sin(t)
#     theta = theta_plot.sin(t)
#     # f = f_fakeValueGenerator.sawtooth(t)
#     # update animation
#     state = np.array([[z], [theta], [0.0], [0.0]])
#     animation.update(state)
#     dataPlot.update(t, state, 1.0) #f)
#     # advance time by t_plot
#     t += P.t_plot
#     plt.pause(0.02)  # allow time for animation to draw

# # Keeps the program from closing until the user presses a button.
# print('Press key to close')
# plt.waitforbuttonpress()
# plt.close()

import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter

#instantiate reference input classes
z_plot = signalGenerator(amplitude=P.length/4., frequency=0.5, y_offset=P.length/2.)
theta_plot = signalGenerator(amplitude=np.pi/8, frequency=0.1, y_offset=0.0)
f_plot = signalGenerator(amplitude=5, frequency=.5)
force = signalGenerator(amplitude=0.5, frequency=1.0, y_offset=11.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
t = P.t_start # time starts at t_start
while t < P.t_end: # main simulation loop
    # set variables
    z = z_plot.sin(t)
    theta = theta_plot.sin(t)
    f = f_plot.sawtooth(t)
    u = force.sin(t)
    # update animation
    state = np.array([[z], [theta], [f_plot], [force]])
    animation.update(state)
    dataPlot.update(t, state, f)
    t += P.t_plot # advance time by t_plot
    plt.pause(0.02)
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()