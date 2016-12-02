import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import csv
import random

#Constants
P = 1				# atm
pi = 3.14159
#R = 4.378 * 10 ** -3	# BTU/(mol*F)
Tref = 77			# F
T_sp = 500			# F

# Grill metrics
V_pot = 0.0818 	# ft^3
V_grill = 5.5 # ft^3
m_steel_fb = 10 #lbs
m_steel_grill = 30 #lbs

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
Cp = np.zeros(8)
Cp[0] = 0.2195		# O2
Cp[1] = 0.2484		# N2
Cp[2] = 0.2016		# CO2
Cp[3] = 0.4705		# H2O
Cp[4] = 0.4777		# Volatiles (guess)	
Cp[5] = 4.2008		# Wood (at 600F)
Cp[6] = 0.2412		# Air
Cp[7] = 0.12

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
#V_air_in = 1		# Change this to correctly be utilized
T_amb = 70			# F		

M = np.zeros(6)
M[0] # O2
M[1] # N2
M[2] # CO2
M[3] # H2O
M[4] # Volatiles	
M[5] # Wood (C32H46O21)	
 

def traeger(X, t, fanSpeed):

    XO2 =  X[0] # Oxygen mass fraction in grill
    XN2 =  X[1] # Nitrogen mass fraction in grill
    XCO2 = X[2] # Carbon Dioxide mass fraction in grill
    XH2O = X[3] # Water Vapor mass fraction in grill
    Temp = X[4]
    T2   = X[5]
    m_wood_in = 1. #lb/hr
    #if i >=100:
     #   m_wood_in = 0.	
    # Rate of combustion of wood --- estimated from what traeger experts say on forums
    rWood = -0.75 * m_wood_in       # [lb/hr]

    # Rate of combustion of O2, H2O, CO2 based off of stoichiometry w/ wood
    rO2 = 1.3786 * rWood         # [lb/hr]
    rCO2 = -1.8381 * rWood        # [lb/hr]
    rH2O = -0.5405 * rWood        # [lb/hr]

    # Air composition
    X_O2_in = 0.21        # Mass fraction of O2 in air
    X_N2_in = 1. - X_O2_in # Mass fraction of N2 in air
    m_grill = V_pot * rho[6] + 1         # Mass of air in grill [lb]

    # Dependent on fan speed
    m_air = 15. * 60. * rho[6] * fanSpeed     #cfm * (60 min/ 1hr) * density * fanSpeed # Flowrate of air (flowrate in and out through flue) [lb/hr]

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
    h_in_air  = (X_O2_in * h_f[0] + X_N2_in * h_f[1]) + (X_O2_in * Cp[0] + X_N2_in * Cp[1]) * (T_amb-Tref)
    h_in_wood = h_f[5] + Cp[5] * (T_amb-Tref)
    h_in_tot  = h_in_air * m_air + h_in_wood * (-rWood)
    
    #Enthalpy of outlet
    h_out = (XO2 * h_f[0] + XN2 * h_f[1] + XCO2 * h_f[2] + XH2O * h_f[3] + Cp_flue * (Temp - Tref)) * m_flue + (m_wood_in + rWood) * Cp[5] *(Temp-Tref)   
    
    #Energy Balance Fire Pot
    dTdt =  (h_in_tot - h_out- 10 *(Temp-70)) * (1 / m_grill) * (1 / Cp_flue)  

    #Energy Balance Grill
    m_grill_tot = m_steel_grill + V_grill * rho[6]
    x_grill_air = (V_grill * rho[6]) / m_grill_tot
    x_grill_steel = m_steel_grill / m_grill_tot
    Cp_grill = Cp[7] * x_grill_steel + Cp[6] * x_grill_air
    dT2dt = (m_flue * Cp_flue *(Temp - T2) - 10 *(T2 - 70)) / Cp_grill / m_grill_tot  
    
    return [dXO2dt, dXN2dt, dXCO2dt, dXH2Odt, dTdt,dT2dt]

# Initial conditions
XO2_0 = 0.2   # Mass fraction
XN2_0 = 0.8   # Mass fractionXCO2_0 = 0.0  # Mass fraction
XH2O_0 = 0.0  # Mass fraction
XCO2_0 = 0.0
Temp_0 = 70  #deg F temp of fire pot
Temp2_0 = 70
X = [XO2_0, XN2_0, XCO2_0, XH2O_0, Temp_0, Temp2_0]

# Time points for solution
tf = 10. # hr
nsteps = 201
t = np.linspace(0, tf, nsteps)
delta_t = tf/(nsteps - 1)

# Fan speed can be adjusted
fan = np.ones((nsteps))
fan[5:50] = 0.5 # Units TBD
fan[50:100] = 1.0
fan[100:150] = 0.1
fan[150:] = 0.5

# Record solution
X_O2        = np.empty((nsteps))
X_O2[0]     = XO2_0
X_N2        = np.empty((nsteps))
X_N2[0]     = XN2_0
X_CO2       = np.empty((nsteps))
X_CO2[0]    = XCO2_0
X_H2O       = np.empty((nsteps))
X_H2O[0]    = XH2O_0
Temp_out    = np.empty((nsteps))
Temp2       = np.empty((nsteps))
Temp2[0]    = Temp2_0
Temp_out[0] = Temp_0

# Simulate the grill step test
for i in range(nsteps-1):
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
    Temp2[i+1]    = result[-1,5]

    # Reset initial conditions
    X[0] = result[-1,0]
    X[1] = result[-1,1]
    X[2] = result[-1,2]
    X[3] = result[-1,3]
    X[4] = result[-1,4]
    X[5] = result[-1,5]
    #print("Sum of mass fractions = ", X[0] + X[1] + X[2] + X[3])
    

# Organize results
data = np.vstack((t,fan,Temp_out,Temp2))
data = data.T
file_name = "data_"+str(int(random.random()*1000))+".txt"
np.savetxt(file_name,data,delimiter=',')
print file_name
# Plot results
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,X_O2[:], 'b-', label="XO2")
plt.plot(t,X_N2[:], 'k:', label="XN2")
plt.plot(t,X_CO2[:], 'k-', label="XCO2")
plt.plot(t,X_H2O[:], 'b:', label="XH2O")
plt.legend(loc="best")
plt.ylabel("O2 Mass Fraction")
plt.xlabel("Time [hr]")

plt.subplot(2,1,2)
plt.plot(t,fan*15*rho[6]*60,'k-')
plt.ylabel("Air Flow Rate [lb/hr]")
plt.xlabel("Time [hr]")
#plt.show()

plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t,Temp_out[:],'b--',label="Temp FireBox")
plt.plot(t,Temp2,'r.',label="Temp Grill")
plt.legend(loc="best")
plt.xlabel("Time")
plt.ylabel("Temp Deg F")
plt.subplot(2,1,2)
plt.plot(t,fan*15*rho[6]*60,'k-')
plt.ylabel("Air Flow Rate [lb/hr]")
plt.xlabel("Time [hr]")
plt.show()
