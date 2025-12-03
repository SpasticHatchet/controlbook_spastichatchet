import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from dataPlotterObserver import dataPlotterObserver
from blockbeamDynamics import blockbeamDynamics
from ctrlLQR import ctrlLQR


# instantiate blockbeam, controller, and reference classes
blockbeam = blockbeamDynamics()
controller = ctrlLQR()
reference = signalGenerator(amplitude=0.125, frequency=0.05, y_offset=0.25)
disturbance = signalGenerator(amplitude=0.25, frequency=0.0)
noise_z = signalGenerator(amplitude=0.01)
noise_th = signalGenerator(amplitude=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = blockbeamAnimation()

t = P.t_start
y = blockbeam.h()
while t < P.t_end:

    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = 0.0 #disturbance.step(t)
        n = np.array([[0.0], [0.0]]) #np.array([[noise_z.random(t)], [noise_th.random(t)]])
        u, xhat = controller.update(r, y + n)
        u = u[0,0]
        y = blockbeam.update(u + d)
        t += P.Ts

    # update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, blockbeam.state, u, r)
    dataPlotObserver.update(t, blockbeam.state, xhat, d, 0.0)
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
