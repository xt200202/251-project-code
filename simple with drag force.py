#!/usr/bin/env python
# coding: utf-8

# In[58]:


import numpy as np
import matplotlib.pyplot as plt
import os

# Constants
g = 9.81  # m/s^2
l = 1.0  # m
dtheta_0_dt = 0.0  # rad/s
m = 1.0  # kg
c = 0.5  # kg/s
tau_drag = 4 * np.pi / np.sqrt((6 * g / l) - (c ** 2) / (m ** 2)) # s oscillation period of the rod pendulum with drag force

# Time settings
t_max = tau_drag*5 # s maximum t 
dt = 0.0001  # s
time = np.arange(0, t_max, dt)

# Initial conditions and damping coefficients
theta_0_values = [0.1, 1.0]
damping_coefficients = [0, 0.5, 1.0]

# Euler method to solve the ODE(Rod pendulum with drag force)
def euler_method(theta_0, damping_coefficient):
    theta = np.zeros(len(time))
    phi = np.zeros(len(time))
    theta[0] = theta_0
    phi[0] = dtheta_0_dt

    for i in range(1, len(time)):
        theta[i] = theta[i - 1] + phi[i - 1] * dt
        phi[i] = phi[i - 1] - (g / (l)) * np.sin(theta[i - 1]) * dt - (damping_coefficient / (m)) * phi[i - 1] * dt

    return theta, phi

# Analytical solution of the problem
def analytical_solution(theta_0):
    return np.exp(alpha * time) * (theta_0 * np.cos(omega * time) - (alpha / omega) * theta_0 * np.sin(omega * time))


# First plot: Comparison of numerical and analytical solutions for different initial conditions
for theta_0 in theta_0_values:
    alpha = -c / (2 * m)
    omega = np.sqrt(g / l - (c / (2 * m)) ** 2)
    theta_numerical, _ = euler_method(theta_0, c)
    theta_analytical = analytical_solution(theta_0)
    plt.plot(time/tau_drag, theta_numerical, label=f"Numerical, θ₀ = {theta_0} rad")
    plt.plot(time/tau_drag, theta_analytical, label=f"Analytical, θ₀ = {theta_0} rad", linestyle="--")

plt.xlabel("Dimensionless Time (t/τ) ")
plt.ylabel("Swing angle (radians)")
plt.legend()
plt.show()

# Second plot: Phase space plots for different damping coefficients
for c in damping_coefficients:
    theta_0 = 1.0
    theta, phi = euler_method(theta_0, c)
    plt.plot(theta, phi, label=f"Damping Coefficient = {c}")

plt.xlabel("Swing Angle (θ) [radians]")
plt.ylabel("Angular Velocity (dθ/dt) [rad/s]")
plt.legend()
plt.grid()
plt.show()


# In[ ]:




