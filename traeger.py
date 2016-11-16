import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

def traeger(temp, t, fanSpeed):
    grillTemp = temp # K

    c1 = 0.5         # Constant that will be changed later for convection

    # Calculate the dTdt from convection
    # dTdt = c1 * fanSpeed * grillTemp
    dTdt = c1 * fanSpeed * (1/grillTemp)

    return dTdt

# Initial conditions
T0 = 300   # K

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
    T = odeint(traeger, T0, [0,1], inputs)
    # Record result
    y[i+1] = T[-1]
    # Reset initial conditions
    T0 = T[-1]

# Plot results
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,y[:], 'b-')
plt.ylabel("Temperature [K]")
plt.xlabel("Time [s]")

plt.subplot(2,1,2)
plt.plot(t,fan,'k-')
plt.ylabel("Fan Speed [ft/s]")
plt.xlabel("Time [s]")
plt.show()
