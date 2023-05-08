#!/usr/bin/env python
# coding: utf-8

# In[19]:


import numpy as np
import matplotlib.pyplot as plt
import os

# Constants
L = 1  # Length of the pendulum in meters
g = 9.81  # Gravitational acceleration in m/s^2
tau = (2*np.pi*np.sqrt((2*L)/(3*g)))  # Time period of the pendulum
theta_0 = 1  # Initial angle in radians
dtheta_0_dt = 0  # Initial angular velocity in radians/s

t_end = 5 * (2*np.pi*np.sqrt((2*L)/(3*g)))  # 5 time periods

# Euler Method for the rod oscillation without drag force
def euler_method_pendulum(dt, initial_angle, t_end):

    theta_vals = []
    theta = initial_angle
    dz = dtheta_0_dt

    t = np.arange(0, t_end, dt)
    for i in range(len(t)):
        theta_vals.append(theta)
        dz_new = dz - ((3*g)/(2*L))*np.sin(theta)*dt
        theta_new = theta + dz * dt
        dz = dz_new
        theta = theta_new

    return t / tau, theta_vals, tau

# The analytical solution of rod pendulum without drag force
def linear_pendulum_analytical(t, initial_angle):
    omega = np.sqrt((3*g)/(2*L))
    return initial_angle * np.cos(omega * t),tau

# plot: Numerical vs Analytical solution
dt = 0.0001
initial_angles = [0.1, 1]

plt.figure()
for angle in initial_angles:
    t_dim, theta_vals, tau = euler_method_pendulum(dt, angle, t_end)
    plt.plot(t_dim, theta_vals, label=f"Numerical, θ₀ = {angle}")

    theta_analytical, tau = linear_pendulum_analytical(t_dim * tau, angle)
    plt.plot(t_dim, theta_analytical, linestyle='--', label=f"Analytical, θ₀ = {angle}")

plt.xlabel("Dimensionless Time ($t/\\tau$)")
plt.ylabel("Swing Angle (radians)")
plt.legend()
plt.show()



# In[ ]:




