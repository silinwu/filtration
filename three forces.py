

#####Extension libraries used in this .py##########################################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import math

R1 = 10e-6
R2 = 10e-6
D = np.arange(1e-13,1e-8,5e-13)
# ELECtrolic force
zeta = 15
fai = 5.4 * zeta
Z = (9.22e-11)*((np.tanh(fai/103))**2)
kappa = 1e9
ELE = kappa*((R1*R2)/(R1+R2))*Z*(np.exp(-kappa*D))
#Van Der fan force
A = 3e-20
VdwCutoff = (A/(12.0*Z*kappa))**0.5
VdwCutoff = float(format(VdwCutoff,'.2g'))
D1 = np.arange(1e-13,VdwCutoff,5e-13)
FWD1 = -(A/(6*VdwCutoff*VdwCutoff))*((R1*R2)/(R1+R2))*D1/D1
D2 = np.arange(VdwCutoff,1e-8,5e-13)
FWD2 = -(A/(6*D2*D2))*((R1*R2)/(R1+R2))
FWD = np.append(FWD1,FWD2)
#Bridging Force
bbeta = 0.1
s = 100e-9
epsilon = 4.11e-20
Lc = 4.28e-6
l = 3.04e-10
BRI = -((Lc-D)/l)*(bbeta*4*math.pi*((R1*R2)/(R1+R2))*epsilon)/(s*s)
#Total force
Force = FWD + BRI + ELE
#Draw
ax = plt.subplot(111)
plt.plot(D,FWD,label="FWD",alpha = 0.6)
plt.plot(D,BRI,label ="BRI",alpha =0.6)
plt.plot(D,ELE,label ="ELE",alpha =0.6)
plt.plot(D,Force,label ="Force",linewidth = 4, alpha =0.6)
plt.legend(loc='upper left', bbox_to_anchor=(0.2,0.6))
ax.set_xscale('log')
plt.hlines(1e-6,1e-13,1e-8)
plt.hlines(0,1e-13,1e-8)
plt.show()


##########################################