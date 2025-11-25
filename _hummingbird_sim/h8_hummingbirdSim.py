import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlPD import ctrlPD

# instantiate pendulum, controller, and reference classes
hummingbird = HummingbirdDynamics(alpha=0.0)
controller = ctrlPD()
disturbance = SignalGenerator(amplitude=0.00)
psi_ref = SignalGenerator(amplitude=np.radians(20), frequency=0.04)
theta_ref = SignalGenerator(amplitude=np.radians(15), frequency=0.05)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = hummingbird.h()
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = np.array([[theta_ref.square(t)], [psi_ref.square(t)]])
        pwm, y_ref = controller.update(r, y)

        d = np.array([  [disturbance.step(t) / P.d],               
                        [disturbance.step(t) / P.d]]) / (2 * P.km)


        y = hummingbird.update(pwm + d)  # Propagate the dynamics
        t += P.Ts  # advance time by Ts

    # update animation and data plots at rate t_plot
    animation.update(t, hummingbird.state)
    dataPlot.update(t, hummingbird.state, pwm, y_ref)

    # the pause causes figure to be displayed during simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
