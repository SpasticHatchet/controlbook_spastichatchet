import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1.5, frequency=0.05)
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.05)
psi_ref = SignalGenerator(amplitude=0.5, frequency=.05)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
dynamics = HummingbirdDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    F = 0.23
    f_l = F/2
    f_r = F/2.2

    u = np.array([[f_l],[f_r]])

    dynamics.update(u)

    state = dynamics.state
    # set variables
    phi = 0 #phi_ref.sin(t) # Rotation of the hummingbird around the main arm.
    theta = 0 # theta_ref.sin(t) # Pitch of the main arm.
    psi = 0 # psi_ref.sin(t)# Rotation of the whole system around z.

    # define dummy state and reference values (since we aren't simulating yet)
    state = np.array([[phi], [theta], [psi], [0.0], [0.0], [0.0]])
    ref = np.array([[0], [0], [0]])
    
    # convert force and torque to pwm values
    force = 0
    torque = 0
    pwm = P.mixing @ np.array([[force], [torque]]) / P.km

    # update animation and data plots
    animation.update(t, state)
    dataPlot.update(t, state, pwm)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()