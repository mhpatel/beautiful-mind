# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 20:04:28 2016

@author: MIT
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import *
'''
In this project, analysis of countercurrent double pipe heat exchanger which used to cool a hot fluid (flow
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
Cph= 1900 # specific heat capacity of oil in J/kg.C
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
L= 4.5 # length of heat exchanger
#correlation for finding heat transfer coefficients:
#for hot fluid:
#Nuh=0.023*Reh**0.8*Prh**0.3
#for Cold fluid:
#Nuc=0.022*Rec**0.8*Prc**0.5
Ac=0.0009425 #m**2
Ah=0.0003141 #m**2
Reh=4*mh/di*0.00001*pi
Rec=4*mc/pi*do*0.001
Prh=140
Prc=4.7
hh=123.64
hc=5421.08
F=hc*Ac/mc*Cpc
G=hh*Ah/mh*Cph
'''N = 100 # number of points to discretize
X = np.linspace(0, L, N) # position along the cylindrical wall
h = L / (N - 1)'''
 
 
#ma=1;mb=2 #ma-mass flow rate of a and mb=mass flow rate of b
P=10.0 #perimeter
U=100.0 #overall heat transfer coefficient
L=10.0 #length of heat exchanger
#n=30.0 #no of points
n=[10,20,30,40,50,60,70,80,90,100]
#p=range(10)
for j in range(10):
    m=0;t1=30;e=0.2;
    while e>0.001:
        to=100;h=0.1;too=t1+m*0.001;
        tB=too;
        for i in range(n[j]):
            tc=too;
            th=to;
            delx=L/(n[j]-1);
            Cph=1900;
            d=-P*U/(mh*Cph);
            dthdx=d*(th-tc);
            to=th+dthdx*delx;
            tA=to;
            Cpc=4200;
            d1=-P*U/(mc*Cpc);
            dtcdx=d1*(th-tc);
            too=tb+dtbdx*delx;
            e=too-30;
            if e<0:
                e=30-too;
                m=m+1;
print th
print tc