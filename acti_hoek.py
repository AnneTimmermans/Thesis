# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:30:46 2020

@author: annee
"""


import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

s = 200
#d = 100
z = np.linspace(100000,1,1000)



hoek = range(70,90,5)
afstanden = range(50, 301, 50)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []

print("L_max,Z_max, acti_L, d")
plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)

    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)
    

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3



    lines0.append(plt.plot(cti/0.3,np.abs(ctpi)/0.3)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines0, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([0,300000])
plt.xlim([0,300])
plt.title("-t' against t for different angles at s=200")
plt.xlabel("t(ns)")
plt.ylabel("-t'(ns)")
plt.show()


plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi

    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3



    #print(L_max,Z_max, acti_L)
    lines1.append(plt.plot(z,np.abs(acti))[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines1, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,80000])
#plt.xlim([0,3])
plt.xlabel("z(m)")
plt.ylabel("Boost")
plt.title("Boost against z with <n(z)> for different angles at s=200")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi

    an = 0 #(k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3



    #print(L_max,Z_max, acti_L)
    lines2.append(plt.plot(z,np.abs(acti))[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines2, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,8000])
#plt.xlim([0,3])
plt.xlabel("z(m)")
plt.ylabel("1/boost")
plt.title("1/boost against z with n(Xmax) for different angles at s=200")
plt.show()

print("L_max,Z_max, acti_L, t met n(z)")
plt.figure(dpi=100)#Gepasseerde massa
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3

    x = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    n_L2 = 1.003 #1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L2 = 0 # (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L2 = an_L2 * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L2 = n_L2 * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L2 = 1
    acti_L2 = np.abs(acti_1_L2 + acti_2_L2 + acti_3_L2)


    print(L_max,Z_max, acti_L, t)
    #lines1.append(plt.plot(cti/0.3,z)[0])
    #lines2.append(plt.plot(acti,z)[0])
    lines3.append(plt.plot(cti/0.3,a)[0])


plt.legend(lines3, hoek)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([0,100000])
plt.xlim([0,100])
plt.title("t against z with <n(z)> for different angles at s=200")
plt.xlabel("t(ns)")
plt.ylabel("z(m)")
plt.show()

print("L_max,Z_max, acti_L, t met n=1.003")
plt.figure(dpi=100)#Gepasseerde massa
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = 0 #(k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3

    x = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    n_L2 = 1.003 #1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L2 = 0 # (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L2 = an_L2 * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L2 = n_L2 * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L2 = 1
    acti_L2 = np.abs(acti_1_L2 + acti_2_L2 + acti_3_L2)


    print(L_max,Z_max, acti_L2, t)
    #lines1.append(plt.plot(cti/0.3,z)[0])
    #lines2.append(plt.plot(acti,z)[0])
    lines4.append(plt.plot(cti/0.3,a)[0])


plt.legend(lines4, hoek)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([0,100000])
plt.xlim([0,100])
plt.title("t against z with n(Xmax) for different angles at s=200")
plt.xlabel("t(ns)")
plt.ylabel("z(m)")
plt.show()