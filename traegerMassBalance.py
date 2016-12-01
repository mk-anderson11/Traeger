import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

def traeger(X, t, fanSpeed):

    XO2 =  X[0] # Oxygen mass fraction in grill
    XN2 =  X[1] # Nitrogen mass fraction in grill
    XCO2 = X[2] # Carbon Dioxide mass fraction in grill
    XH2O = X[3] # Water Vapor mass fraction in grill


    # Rate of combustion of wood --- estimated from what traeger experts say on forums
    rWood = -0.5       # [lb/hr]

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
    m_flue = m_air + (-rWood)      # Flowrate of flue vapor [lb/hr]

    # Calculate the dXO2dt
    # dXO2dt = (1/m_grill) * (m_air * (X_O2_in - XO2) + rO2)
    dXO2dt  = (1/m_grill)*(X_O2_in*m_air - XO2*m_flue + rO2)
    dXN2dt  = (1/m_grill)*(X_N2_in*m_air - XN2*m_flue)
    dXCO2dt = (1/m_grill)*(-XCO2*m_flue + rCO2)
    dXH2Odt = (1/m_grill)*(-XH2O*m_flue + rH2O)

    return [dXO2dt, dXN2dt, dXCO2dt, dXH2Odt]

# Initial conditions
XO2_0 = 0.3   # Mass fraction
XN2_0 = 0.6   # Mass fraction
XCO2_0 = 0.1  # Mass fraction
XH2O_0 = 0.0  # Mass fraction
X = [XO2_0, XN2_0, XCO2_0, XH2O_0]

# Time points for solution
tf = 10 # hr
nsteps = 11
t = np.linspace(0, tf, tf+1)
delta_t = tf/(nsteps - 1)

# Fan speed can be adjusted
fan = np.ones((tf+1))
fan[5:50] = 5 # Units TBD
fan[50:100] = 10
fan[100:150] = 2
fan[150:] = 5

# Record solution
X_O2     = np.empty((tf+1))
X_O2[0]  = XO2_0
X_N2     = np.empty((tf+1))
X_N2[0]  = XN2_0
X_CO2    = np.empty((tf+1))
X_CO2[0] = XCO2_0
X_H2O    = np.empty((tf+1))
X_H2O[0] = XH2O_0

# Simulate the grill step test
for i in range(tf):
    # Specify fan speed
    inputs = (fan[i],)
    # Integrate the model

    # XO2 = odeint(traeger, XO2_0, [0,delta_t], inputs)
    result = odeint(traeger, X, [0,delta_t], inputs)

    # Record result
    X_O2[i+1]  = result[-1,0]
    X_N2[i+1]  = result[-1,1]
    X_CO2[i+1] = result[-1,2]
    X_H2O[i+1] = result[-1,3]

    # Reset initial conditions
    X[0] = result[-1,0]
    X[1] = result[-1,1]
    X[2] = result[-1,2]
    X[3] = result[-1,3]

    print("Sum of mass fractions = ", X[0] + X[1] + X[2] + X[3])

# Organize results


# Plot results
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,X_O2[:], 'b-', label="XO2")
plt.plot(t,X_N2[:], 'k:', label="XN2")
plt.plot(t,X_CO2[:], 'k-', label="XCO2")
plt.plot(t,X_H2O[:], 'b:', label="XH2O")
plt.legend(loc="best")
plt.ylabel("O2 Mass Fraction")
plt.xlabel("Time [s]")

plt.subplot(2,1,2)
plt.plot(t,fan,'k-')
plt.ylabel("Fan Speed [ft/s]")
plt.xlabel("Time [s]")
plt.show()
