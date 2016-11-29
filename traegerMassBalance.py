import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

def traeger(X, t, fanSpeed):

    XO2 = X # Oxygen mass in grill [lb]



    # Rate of combustion of wood --- estimated from what traeger experts say on forums
    rWood = -1.0       # [lb/hr]

    # Rate of combustion of O2, H2O, CO2 based off of stoichiometry w/ wood
    rO2 = 1.3786 * rWood         # [lb/hr]
    rCO2 = -1.8381 * rWood        # [lb/hr]
    rH2O = -0.5405 * rWood        # [lb/hr]

    # Air composition
    X_O2_in = 0.21        # Mass fraction of O2 in air
    X_N2_in = 1 - X_O2_in # Mass fraction of N2 in air
    m_grill = 10.         # Mass of air in grill [lb]

    # Dependent on fan speed
    m_air = 8.*fanSpeed     # Flowrate of air (flowrate in and out through flue) [lb/hr]

    # Dependent on combustion
    m_flue = 5.             # Flowrate of flue vapor [lb/hr]

    # Calculate the dXO2dt
    dXO2dt = (1/m_grill) * (m_air * (X_O2_in - XO2) + rO2)

    return dXO2dt

# Initial conditions
XO2_0 = 0.5   # K

# Time points for solution
tf = 10 # hr
nsteps = 11
t = np.linspace(0, tf, tf+1)
delta_t = tf/(nsteps - 1)

# Fan speed can be adjusted
fan = np.ones((tf+1))
fan[10:50] = 5 # Units TBD
fan[50:100] = 10
fan[100:150] = 2
fan[150:] = 5

# Record solution
y = np.empty((tf+1))
y[0] = XO2_0

# Simulate the grill step test
for i in range(tf):
    # Specify fan speed
    inputs = (fan[i],)
    # Integrate the model

    print("XO2_0", XO2_0)

    XO2 = odeint(traeger, XO2_0, [0,delta_t], inputs)

    print("XO2 = ", XO2[-1][0])

    # Record result
    y[i+1] = XO2[-1]
    # Reset initial conditions
    XO2_0 = XO2[-1]

# Plot results
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,y[:], 'b-')
plt.ylabel("O2 Mass Fraction")
plt.xlabel("Time [s]")

plt.subplot(2,1,2)
plt.plot(t,fan,'k-')
plt.ylabel("Fan Speed [ft/s]")
plt.xlabel("Time [s]")
plt.show()
