# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 10:21:59 2016

@author: user
"""
from scipy.integrate import odeint
from scipy import linspace
from scipy import array
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy.optimize as optimization
kla=0.07
A=10
pT=1#"""atm"""
Hc=29.4#"""atm/M"""
MDc=4.089*10**(-2)#"M"
Hm=714.286#"atm/M"
MDm=4.1*10**(-2)#"M"
kGa=7*10**(-7)
Pwsat=0.03125#"atm"
def derivate(y,t):
    yc=y[0]/(y[0]+y[1]+y[2])
    ym=y[1]/(y[0]+y[1]+y[2])
    yw=y[2]/(y[0]+y[1]+y[2])
    xc=y[3]/(y[3]+y[4]+y[5])
    xm=y[4]/(y[3]+y[4]+y[5])
    #xw=y[5]/(y[3]+y[4]+y[5])
    dGc=-kla*A*(yc*pT/Hc-MDc*xc)
    dLc=-kla*A*(yc*pT/Hc-MDc*xc)
    dGm=-kla*A*(ym*pT/Hm-MDm*xm)
    dLm=-kla*A*(ym*pT/Hm-MDm*xm)
    dGw=kGa*A*(pT*yw-Pwsat)
    dLw=kGa*A*(pT*yw-Pwsat)
    return array([dGc,dLc,dGm,dLm,dGw,dLw])
a=array ([10,10,80])
t=linspace(0.0,10.0,100)
def error(x,a,b,c):
    #t=linspace(0.0,10.0,100)
    yinitial=([50,a,50,b,0,c])
    y=odeint(derivate,yinitial,t)
    return array([y[99,5],y[99,3],y[99,1]])
#ans=fsolve(error,a)
xdata=array([10])
ydat=array([100,0,0])
ydata=np.transpose(ydat)
print ydata
x0=array([10,10,80])
#sigma=
([ans,err])=optimization.curve_fit(error,xdata,ydata,x0)
print ans