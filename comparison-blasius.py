#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.interpolate import interp1d

#===============================================================
# Defining derivative function
#===============================================================

def deriv(y, eta):
    return np.array([y[1], y[2], -y[0]*y[2]])

#===============================================================
# Defining function to solve ODE
#===============================================================

def solve_ode(fppp0):
    eta = np.linspace(0.0, 6.0, 1000)
    y0 = np.array([0.0, 0.0, fppp0])
    y = odeint(deriv, y0, eta)
    return (y[-1,1]-1.0)**2

#===============================================================
# Setting defaults
#===============================================================

lwidth = 1.5

#===============================================================
# Inputs
#===============================================================

nu = 15.61434409e-6     # m^2/s
Ue = 5.0                # m/s
x = 0.3                 # m

# Running "shooting" procedure to determine which boundary 
# condition will yield the correct profile
#res = minimize(solve_ode, 2.0)
res = minimize(solve_ode, 1.25)
print("res = {}".format(res))

# Solving the ODE again with the correct boundary condition
eta = np.linspace(0.0, 4.5, 4501)
y0 = np.array([0.0, 0.0, res.x])
y = odeint(deriv, y0, eta)

# Importing numerical data from file
d = np.genfromtxt("explicit_profile.dat", skip_header=2)
#d = np.genfromtxt("crank-nicolson-profile.dat", skip_header=2)
y_n = d[:, 0]
u_n = d[:, 1]
v_n = d[:, 2]

# Converting to eta
eta_n = y_n*np.sqrt(Ue/(2.0*nu*x))

# Separating values from the solution
f = y[:,0]
fp = y[:,1]
fpp = y[:,2]

#===============================================================
# Plotting
#===============================================================

# Plotting velocity profiles
fig1 = plt.figure(figsize=(6.5,4.5))
plt.plot(eta, fp, "-k", linewidth=lwidth, label="Blasius")
plt.plot(eta_n, u_n/Ue, "--r", lw=lwidth, label="Numerical")
plt.xlabel(r"$\eta$")
plt.ylabel(r"$u/U_e$")
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 4.5])
plt.legend(loc=2)
plt.tight_layout()

# Making sure images directory exists
if not os.path.isdir("images"):
    os.mkdir("images")

# Saving figure
fig1.savefig("images/u.pdf")

# Determining v(y) at x = 0.3 m
Re_x = Ue*x/nu
delta = 5.0*x/np.sqrt(Re_x)
#n = np.linspace(0.0, delta, 500)
#y = odeint(deriv, y0, n)
#fn = y[:,0] 
#fpn = y[:,1]
v = np.sqrt(nu*Ue/(2.0*x))*(eta*fp - f)
fig2 = plt.figure(figsize=(6,4))
plt.plot(v/Ue, eta, "-k", linewidth=lwidth, label="Blasius")
plt.plot(v_n/Ue, eta_n, "--r", linewidth=lwidth, label="Numerical")
plt.xlabel(r"$v/U_e$")
plt.ylabel(r"$\eta$")
plt.ylim([0.0, 4.5])
plt.legend(loc=2)
plt.tight_layout()

# Showing figure
fig2.savefig("images/v.pdf")
plt.show()


