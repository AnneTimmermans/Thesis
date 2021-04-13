# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 20:11:53 2020

@author: annee
"""


import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#s = 100
#d = 100
z = np.linspace(100000,1,10000)



hoek = range(0,90,15)
afstanden = range(100, 1000, 200)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []


plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
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
    D = acti
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    #Number of particles distribution N
    Ep = 10**19 # eV
    Ne = 6 * (Ep / 10**10)
    
    X_0 = 36.7# * 10**4 # g/m2
    S = (3*X/X_0)/(X/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    A_0 = (Cx*N)/D
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = A_0/A_x#Cx/v


    print(Ft)
    lines0.append(plt.plot(a,Cx)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend(lines0, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,2000])
#plt.xlim([0,400])
plt.title("Cx against z")
plt.xlabel("z(m)")
plt.ylabel("Cx")
plt.show()

plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
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
    D = acti
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    #Number of particles distribution N
    Ep = 10**19 # eV
    Ne = 6 * (Ep / 10**10)
    
    X_0 = 36.7# * 10**4 # g/m2
    S = (3*X/X_0)/(X/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    A_0 = (Cx*N)/D
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v



    lines1.append(plt.plot(a,N)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend(lines1, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,1e61])
plt.xlim([0,50000])
plt.title("N against z")
plt.xlabel("z(m)")
plt.ylabel("N")
plt.show()

plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
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
    D = acti
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    #Number of particles distribution N
    Ep = 10**19 # eV
    Ne = 6 * (Ep / 10**10)
    
    X_0 = 36.7# * 10**4 # g/m2
    S = (3*X/X_0)/(X/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    A_0 = (Cx*N)/D
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v



    lines2.append(plt.plot((cti)/0.3,D)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend(lines2, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,1.2e44])
plt.xlim([0,400])
plt.title("D against t")
plt.xlabel("t(ns)")
plt.ylabel("D")
plt.show()

plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168 # g m-3
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
    D = acti
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    #Number of particles distribution N
    Ep = 10**19 # eV
    Ne = 6 * (Ep / 10**10)
    
    X_0 = 36.7# * 10**4 # g/m2
    S = (3*X/X_0)/(X/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    A_0 = (Cx*N)/D
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v



    lines3.append(plt.plot(rho,a)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend(lines3, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
#plt.xlim([0,1.003])
#plt.ylim([0,10000])
plt.title("rho against z")
plt.xlabel("rho")
plt.ylabel("z")
plt.show()

plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168 # g m-3
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #g/m-3
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    D = acti
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    #Number of particles distribution N
    Ep = 10**19 # eV
    Ne = 6 * (Ep / 10**10)
    
    X_0 = 36.7# * 10**4 # g/m2
    S = (3*X/X_0)/(X/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    A_0 = (Cx*N)/D
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v



    lines3.append(plt.plot(a,v)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend(lines3, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,30000])
#plt.ylim([0,10000])
plt.title("v against z")
plt.xlabel("z(m)")
plt.ylabel("v")
plt.show()

plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168 # g m-3
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #g/m-3
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    D = acti
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    n_L = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    an_L = (k/((-(-L_max))**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_L = an_L * ((-(-L_max))**2 + s**2)**0.5
    acti_2_L = n_L * ((-L_max) / (((-(-L_max))**2 + s**2)**0.5))
    acti_3_L = 1
    acti_L = np.abs(acti_1_L + acti_2_L + acti_3_L)
    
    #Number of particles distribution N
    Ep = 10**19 # eV
    Ne = 6 * (Ep / 10**10)
    
    X_0 = 36.7# * 10**4 # g/m2
    S = (3*X/X_0)/(X/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    A_0 = (Cx*N)/D
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v



    lines4.append(plt.plot((ctpi)/300,v)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend(lines4, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.xlim([-45,0])
#plt.ylim([0,10000])
plt.title("v against t")
plt.xlabel("t(ns)")
plt.ylabel("v")
plt.show()


