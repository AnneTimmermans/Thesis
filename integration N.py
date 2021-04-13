# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:43:46 2020

@author: annee
"""

import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad

s = 200
z = np.linspace(100000,1,1000)



hoek = range(0,90,15)
afstanden = range(100, 2200, 400)
lines1 = []
lines2 = []
lines3 = []
lines4 = []
N = []

def n1(l,R,theta,rho_0,C):
    return 1 + 0.226 * rho_0 * np.exp(-C*(-R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5)) *10**(-6)

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    N = []
    N_2 = []
    delta = 1
    d = s/(np.cos(theta))
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5
#    rho = rho_0 * np.exp(-C*a) #meter
#    n = 1 + 0.226 * rho *10**(-6) # 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    for i in l:
        N.append(quad(n1, i, i+delta, args=(R,theta,rho_0,C))[0])
        N_2.append(quad(n1, i, i+delta, args=(R,theta,rho_0,C))[0]+0.00000000000000001)
    
    n_1 = N
    n_2 = N_2
    ctpi_1 = -l
    ctpi_2 = -l + delta

    cti_1 = n_1 * ((-ctpi_1)**2 + s**2)**0.5 + ctpi_1
    cti_2 = n_2 * ((-ctpi_2)**2 + s**2)**0.5 + ctpi_2
    
    acti_num = (cti_1 - cti_2)/delta
    
    lines1.append(plt.plot(a,N)[0])
    
plt.legend(lines1, afstanden)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,50000])
plt.ylim([1,1.0004])
#plt.xlim([20,140])
#plt.ylim([0.99999,1.00001])
plt.xlabel("z(m)")
plt.ylabel("n")
plt.title(" <n> numerically")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    delta = 1
    d = s/(np.cos(theta))
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5
    k = 0.226*rho_0*10**(-6)
#    rho = rho_0 * np.exp(-C*a) #meter
#    n = 1 + 0.226 * rho *10**(-6) # 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti_ana = acti_1 + acti_2 + acti_3
    
    
    lines2.append(plt.plot(a,n)[0])
    
plt.legend(lines2, afstanden)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,50000])
#plt.xlim([20,140])
plt.ylim([1,1.0004])
#plt.ylim([0.99999,1.00001])
plt.xlabel("z(m)")
plt.ylabel("n")
plt.title("<n> analytically")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    N = []
    N_2 = []
    delta = 1
    d = s/(np.cos(theta))
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5
#    rho = rho_0 * np.exp(-C*a) #meter
#    n = 1 + 0.226 * rho *10**(-6) # 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    for i in l:
        N.append(quad(n1, i, i+delta, args=(R,theta,rho_0,C))[0])
        N_2.append(quad(n1, i, i+delta, args=(R,theta,rho_0,C))[0]+0.00000000000000001)
    
    n_1 = N
    n_2 = N_2
    ctpi_1 = -l
    ctpi_2 = -l + delta

    cti_1 = n_1 * ((-ctpi_1)**2 + s**2)**0.5 + ctpi_1
    cti_2 = n_2 * ((-ctpi_2)**2 + s**2)**0.5 + ctpi_2
    
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    
    acti_ana = acti_1 + acti_2 + acti_3
    acti_num = (cti_1 - cti_2)/delta
    
    acti = acti_num/acti_ana
    acti_diff = acti_num-acti_ana
    
    n_ratio = (N/n)
    n_diff = N-n
    
    lines1.append(plt.plot(a,np.abs(n_diff))[0])
    
plt.legend(lines1, afstanden)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,50000])
#plt.ylim([0.9999,1.0001])
#plt.xlim([20,140])
plt.ylim([0,0.0001])
#plt.ylim([1e-10,1e0])
plt.xlabel("z(m)")
plt.ylabel("Difference")
plt.title("Difference <n> numerical - analytical")
plt.show()