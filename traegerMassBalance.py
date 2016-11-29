import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

def traeger(X, t, fanSpeed):
    XO2 = X # Oxygen mass in grill [lb]

    rO2 = 0.1         # Rate of reaction/combustion O2 [lb/hr]
    m_air = 5         # Flowrate of air (flowrate in and out through flue) [lb/hr]
    X_O2_in = 0.21    # Mass fraction of O2 in air
    m_grill = 10      # Mass of air in grill [lb]

    # Calculate the dXO2dt
    dXO2dt = (1/m_grill) * (m_air * (X_O2_in - XO2) - rO2)

    return dXO2dt

# Initial conditions
XO2_0 = 0.5   # K

# Time points for solution
tf = 200 # s
t = np.linspace(0, tf, tf+1)

# Fan speed can be adjusted
fan = np.ones((tf+1))
fan[10:50] = 5 # Units TBD
fan[50:100] = 10
fan[100:150] = 2
fan[150:] = 5

# Record solution
y = np.empty((tf+1))
y[0] = T0

# Simulate the grill step test
for i in range(tf):
    # Specify fan speed
    inputs = (fan[i],)
    # Integrate the model
    XO2 = odeint(traeger, T0, [0,1], inputs)
    # Record result
    y[i+1] = XO2[-1]
    # Reset initial conditions
    T0 = XO2[-1]

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
