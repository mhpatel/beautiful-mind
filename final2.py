# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 15:09:26 2016

@author:Mit
"""

import numpy as np
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import pylab as plt
from scipy import *
'''
In this project, analysis of co-current double pipe heat exchanger which used to cool a hot fluid (flow
inside tube) by a cooled fluid (flow inside the annulus area) is done. The finite difference techniques was
used to solve the unsteady state countercurrent double pipe heat exchanger differential equations.
The assumptions that used in the theoretical analysis of countercurrent double pipe heat
exchanger are following:
1- There is no phase change for the hot & cold fluids during the heat exchanger pipes.
2- Hot & cold fluids are incompressible & turbulent in flow.
3- The axial heat conduction inside the tubes wall is negligible ( the convection resistance is
   too high compared width conduction resistance).
4- The hot, cold fluids & the pipe wall temperature subject to unsteady state behavior.
5- The specific heat coefficient is indepedent of temperature.
6- The hot fluid flow inside the inner pipe & the cold fluid flow inside the annulus area between pipes.
'''
# Information of Heat Exchanger Fluids:
mh= 0.5 # mass flow rate of hot fluid (oil) in kg/s
mc= 0.3 # mass flow rate of cold fluid (water) in kg/s
Cph= 6900 # specific heat capacity of oil in J/kg.C
Cpc= 4200 # specific heat capacity of water in J/kg.C
vh= 10**-5 # kinematic vicosity of oil in m**2/s
vc= 7*10**-7 # kinematic viscosity of water in m**2/s
kh= 0.134 # thermal conductivity of oil in W/m.C
kc= 0.64 # thermal conductivity of water in W/m.C
Tih= 100 # inlet temperature of oil in degree C
Tic= 30 # inlet temperature of water in degree C

#Information of DOUBLE PIPE Heat Exchanger:
kwi= 384 # thermal conductivity of inner tube
kwo= 45 # thermal conductivity of outer tube
di= 0.02 # inner tube diameter
do= 0.04 # outer tube diameter
ti= 10**-3 # inner tube thickness
to= 3*10**-3 # outer tube thickness
#L= 4.5 # length of heat exchanger
#correlation for finding heat transfer coefficients:
#for hot fluid:
#Nuh=0.023*Reh**0.8*Prh**0.3
#for Cold fluid:
#Nuc=0.022*Rec**0.8*Prc**0.5
Ac=0.0009425 #m**2 heat transfer area for cold fluid
Ah=0.0003141 #m**2 heat transfer area for hot fluid
Reh=4*mh/di*0.00001*pi
Rec=4*mc/pi*do*0.001
Prh=140
Prc=4.7
hh=123.64#W/m2k heat transfer coefficient of hot fluid
hc=5421.08#W/m2k heat transfer coefficient of cold fluid
F=hc*Ac/mc*Cpc
G=hh*Ah/mh*Cph
'''N = 100 # number of points to discretize
X = np.linspace(0, L, N) # position along the cylindrical wall
h = L / (N - 1)'''
 
 
#ma=1;mb=2 #ma-mass flow rate of a and mb=mass flow rate of b
P=10.0 #perimeter
U=10.456 #overall heat transfer coefficient
L=10.0 #length of heat exchanger
dx=10.0/10
tw=np.zeros((11,101))
tc=np.zeros((11,101))
th=np.zeros((11,101))
dtwdx=np.zeros((11,101))
dthdx=np.zeros((11,101))
dtcdx=np.zeros((11,101))
dtwdt=np.zeros((11,101))
dtcdt=np.zeros((11,101))
dthdt=np.zeros((11,101))
th[0,0]=100
e=10
tc[0,0]=40
k1=0.8750
k2=0.1250
F=0.0185#constants
G=0.0018493#constants
a=0
'''
i is the space variable
j is the time variable
'''
for j in range(100):
    
    
    
    for i in range(10):
        dtcdx[i,j]=(P*U*(th[i,j]-tc[i,j]))/(mc*Cpc)
        dthdx[i,j]=(P*U*(th[i,j]-tc[i,j]))/(mh*Cph)
        th[i+1,j]=th[i,j]-dthdx[i,j]*dx
        tc[i+1,j]=tc[i,0]+dtcdx[i,j]*dx
        tw[i,j]=k1*th[i,j]-k2*25.0
        dtcdt[i,j]=(F*(tw[i,j]-tc[i,j])-G*dtcdx[i,j])*scipy.exp(-j*k2*k1)
        tc[i,j+1]=tc[i,j]+1*dtcdt[i,j]
        dthdt[i,j]=(F*(tw[i,j]-th[i,j])+G*dthdx[i,j])*scipy.exp(-j*k1*k2)
        th[i,j+1]=th[i,j]+1*dthdt[i,j]
        dtwdt[i,j]=(F*(th[i,j]-tw[i,j])+G*(tc[i,j]-tw[i,j]))*scipy.exp(-j*k1*k2)
        tw[i,j+1]=tw[i,j]-1*dtwdt[i,j]
    plt.plot(j,th[10,j],'ro')
    #plt.plot(j,tc[10,j])
    #fig=plt.figure()
    #ax=fig.add_subplot(111)
    #plt.plot(j,th[10,j],'ro')
    #ax.title.set_text('temp Vs time')
    #ax.xaxis.label.set_text('time')
    #ax.yaxis.label.set_text('temp')
plt.show()

print 'inlet temperature of cold ',tc[0,0]
print 'inlet temperature of hot ',th[0,0]
print 'outlet temperature of hot ',th[10,99]
print 'outlet temperature of cold ',tc[10,99]