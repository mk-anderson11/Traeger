import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv

#Constants
P = 1				# atm
pi = 3.14159
#R = 4.378 * 10 ** -3	# BTU/(mol*F)
Tref = 100			# F
T_sp = 500			# F

# Grill metrics
V_pot = 0.0818		# ft^3


# Molecular weights (lb/mol)
mw = np.zeros(7)
mw[0] = 0.0705		# O2
mw[1] = 0.0617		# N2
mw[2] = 0.0970		# CO2
mw[3] = 0.0397		# H2O
mw[4] = 1.6887		# Volatiles
mw[5] = 1.6887		# Wood (C32H46O21)
mw[6] = 0.0636		# Air

# Heats of formation (BTU/lb)
h_f = np.zeros(7)
h_f[0] = 0
h_f[1] = 0
h_f[2] = -3844.964
h_f[3] = -5775.724
h_f[4] = 2278.0
h_f[5] = -2278.0
h_f[6] = 0

# Constant heat capacities (BTU/(lb*F))
Cp = np.zeros(7)
Cp[0] = 0.2195		# O2
Cp[1] = 0.2484		# N2
Cp[2] = 0.2016		# CO2
Cp[3] = 0.4705		# H2O
Cp[4] = 0.4777		# Volatiles (guess)
Cp[5] = 4.2008		# Wood (at 600F)
Cp[6] = 0.2412		# Air

# Constant densities (lb/ft^3)
rho = np.zeros(7)
rho[0] = 0.0892		# O2
rho[1] = 0.0727		# N2
rho[2] = 0.1150		# CO2
rho[3] = 0.0480		# H2O
rho[4] = 0.4777		# Volatiles (guess)
rho[5] = 70.757		# Wood
rho[6] = 0.0752		# Air

# Stoichiometry
x = np.zeros(5)
x[0] = -33			# O2    ****(Was originally -31)
x[1] = 0			# N2
x[2] = 32			# CO2
x[3] = 23			# H2O
x[4] = -1			# Volatiles		(added a negative sign)

# Moved from parameters
V_air_in = 1		# Change this to correctly be utilized
T_amb = 70			# F

M = np.zeros(6)
M[0] # O2
M[1] # N2
M[2] # CO2
M[3] # H2O
M[4] # Volatiles
M[5] # Wood (C32H46O21)

counter = 0

def traeger(X, t, fanSpeed):

	XO2 =  X[0] # Oxygen mass fraction in grill
	XN2 =  X[1] # Nitrogen mass fraction in grill
	XCO2 = X[2] # Carbon Dioxide mass fraction in grill
	XH2O = X[3] # Water Vapor mass fraction in grill
	Temp = X[4]

    # Rate of combustion of wood --- estimated from what traeger experts say on forums
	rWood = -1/200       # [lb/hr]

    # Rate of combustion of O2, H2O, CO2 based off of stoichiometry w/ wood
	rO2 = 1.3786 * rWood         # [lb/hr]
	rCO2 = -1.8381 * rWood        # [lb/hr]
	rH2O = -0.5405 * rWood        # [lb/hr]

    # Air composition
	X_O2_in = 0.21        # Mass fraction of O2 in air
	X_N2_in = 1 - X_O2_in # Mass fraction of N2 in air
	m_grill = 1.         # Mass of air in grill [lb]

    # Dependent on fan speed
	m_air = 20 * fanSpeed     # Flowrate of air (flowrate in and out through flue) [lb/hr]

    # Dependent on combustion
	m_flue = m_air + (-rWood)      # Flowrate of flue vapor [lb/hr]

    # Calculate the dXO2dt
    # dXO2dt = (1/m_grill) * (m_air * (X_O2_in - XO2) + rO2)
	dXO2dt  = (1/m_grill)*(X_O2_in*m_air - XO2*m_flue + rO2)
	dXN2dt  = (1/m_grill)*(X_N2_in*m_air - XN2*m_flue)
	dXCO2dt = (1/m_grill)*(-XCO2*m_flue + rCO2)
	dXH2Odt = (1/m_grill)*(-XH2O*m_flue + rH2O)

    #Intermediate Values for Energy Balance

	Cp_flue = XO2 * Cp[0] + XN2 * Cp[1] + XCO2 * Cp[2] + XH2O * Cp[3]  #Heat Capacity of flue gas

    #Enthalpy of inlet (weighted with mass fractions)
	h_in_air  = (X_O2_in * h_f[1] + X_N2_in * h_f[2]) + (X_O2_in * Cp[1] + X_N2_in * Cp[2]) * (T_amb-Tref)
	h_in_wood = h_f[5] + Cp[5] * (T_amb-Tref)
	h_in_tot  = h_in_air * m_air + h_in_wood * (-rWood)

    #Enthalpy of outlet
	h_out = (XO2 * h_f[0] + XN2 * h_f[1] + XCO2 * h_f[2] + XH2O * h_f[3] + Cp_flue * (Temp - Tref)) * m_flue
    #Energy Balance
	dTdt = (h_in_tot - h_out) * (1 / m_grill) * (1 / Cp_flue)
    #print Cp_flue

	# Error checking

	print("h_in_air = ", h_in_air)
	print("dTdt = ", dTdt)


	return [dXO2dt, dXN2dt, dXCO2dt, dXH2Odt, dTdt]

# Initial conditions
XO2_0 = 0.3   # Mass fraction
XN2_0 = 0.6   # Mass fraction
XCO2_0 = 0.1  # Mass fraction
XH2O_0 = 0.0  # Mass fraction
Temp_0 = 200  #deg F temp of fire pot
X = [XO2_0, XN2_0, XCO2_0, XH2O_0, Temp_0]

# Time points for solution
tf = 200 # hr
nsteps = 200
t = np.linspace(0, tf, tf+1)
delta_t = tf/(nsteps - 1)

# Fan speed can be adjusted
fan = np.ones((tf+1))
fan[5:50] = 0.5 # Units TBD
fan[50:100] = 1.0
fan[100:150] = 0.4
fan[150:] = 0.5

# Record solution
X_O2        = np.empty((tf+1))
X_O2[0]     = XO2_0
X_N2        = np.empty((tf+1))
X_N2[0]     = XN2_0
X_CO2       = np.empty((tf+1))
X_CO2[0]    = XCO2_0
X_H2O       = np.empty((tf+1))
X_H2O[0]    = XH2O_0
Temp_out    = np.empty((tf+1))
Temp_out[0] = Temp_0

# Simulate the grill step test
for i in range(tf):
    # Specify fan speed
    inputs = (fan[i],)
    # Integrate the model

    # XO2 = odeint(traeger, XO2_0, [0,delta_t], inputs)
    result = odeint(traeger, X, [0,delta_t], inputs)

    # Record result
    X_O2[i+1]     = result[-1,0]
    X_N2[i+1]     = result[-1,1]
    X_CO2[i+1]    = result[-1,2]
    X_H2O[i+1]    = result[-1,3]
    Temp_out[i+1] = result[-1,4]

    # Reset initial conditions
    X[0] = result[-1,0]
    X[1] = result[-1,1]
    X[2] = result[-1,2]
    X[3] = result[-1,3]
    X[4] = result[-1,4]
    #print("Sum of mass fractions = ", X[0] + X[1] + X[2] + X[3])

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

plt.figure(1)
plt.plot(t,Temp_out[:],'b--',label="Temp")
plt.legend(loc="best")
plt.xlabel("Time")
plt.ylabel("Temp Deg F")
plt.show()
