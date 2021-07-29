# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:53:40 2021

@author: annee
"""

import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#s = 100
#d = 100
z = np.linspace(100000,1,1000)

c = 3*10**8


hoek = range(0,90,15)
afstanden = [1400]#range(100, 2000, 400)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []
lines5 = []


plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 85 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_Xmax = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (L_max**2 + R**2 + 2*R*L_max* np.cos(theta))**0.5)) *10**(-6)
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    
    boost_n = 1 - n * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    boost_nXmax = 1 - n_Xmax * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    boost_ngem = 1 - n_gem * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    
    R_path = (cti-ctpi)/n_gem
    
    #Boost potential
    D_1 = R_path * acti
    D_n = R_path * boost_n
    D_nXmax = R_path * boost_nXmax
    D_ngem = R_path * boost_ngem
    
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
    A_01 = (Cx*N)/D_1
    A_0n = (Cx*N)/D_n
    A_0nXmax = (Cx*N)/D_nXmax
    A_0ngem = (Cx*N)/D_ngem
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    q = 1.6 * 10**(-19)
    schaling = q/(4*np.pi*e_0)
    A_01 = schaling * A_01
    A_0n = schaling * A_0n
    A_0nXmax = schaling * A_0nXmax
    A_0ngem = schaling * A_0ngem
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04 #c
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x1 = (v*N)/D_1
    A_xn = (v*N)/D_n
    A_xnXmax = (v*N)/D_nXmax
    A_xngem = (v*N)/D_ngem
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    q = 1.6 * 10**(-19)
    schaling_x = q/(4*np.pi*e_0*c**2)
    A_x1 = schaling_x * A_x1
    A_xn = schaling_x * A_xn
    A_xnXmax = schaling_x * A_xnXmax
    A_xngem = schaling_x * A_xngem




    lines0.append(plt.plot(a,np.abs(A_01), label= 'dt/dtr')[0])
    lines1.append(plt.plot(a,np.abs(A_0n), label= '1 - n cos(theta)')[0])
    #lines2.append(plt.plot(a,np.abs(A_0ngem), label= '1 - <n> cos(theta)')[0])

    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend()
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([0, 2e-13])
plt.xlim([20000,42000])
plt.title("Charge excess potential against z, s = 1400, theta=85 ")
plt.xlabel("z(m)")
plt.ylabel("A0")
plt.show()

plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 85 * (np.pi/180)
    
    h = 1 # m --> hoogte pannenkoek
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    X_max = 700 #* 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_Xmax = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (L_max**2 + R**2 + 2*R*L_max* np.cos(theta))**0.5)) *10**(-6)
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    
    boost_n = 1 - n*(-ctpi)/(np.sqrt(ctpi**2 + s**2))
    boost_nXmax = 1 - n_Xmax * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    boost_ngem = 1 - n_gem * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    
    R_path = (cti-ctpi)/n_gem
    
    #Boost potential
    D_1 = R_path * acti
    D_n = R_path * boost_n
    D_nXmax = R_path * boost_nXmax
    D_ngem = R_path * boost_ngem
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter


    #Bepaling shower maximum L
    X_max = 700 #* 10**4 g/m2
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
    A_01 = (Cx*N)/D_1
    A_0n = (Cx*N)/D_n
    A_0nXmax = (Cx*N)/D_nXmax
    A_0ngem = (Cx*N)/D_ngem
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    q = 1.6 * 10**(-19)
    schaling = q/(4*np.pi*e_0)
    A_01 = schaling * A_01
    A_0n = schaling * A_0n
    A_0nXmax = schaling * A_0nXmax
    A_0ngem = schaling * A_0ngem
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x1 = (v*N)/D_1
    A_xn = (v*N)/D_n
    A_xnXmax = (v*N)/D_nXmax
    A_xngem = (v*N)/D_ngem
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    q = 1.6 * 10**(-19)
    schaling_x = q/(4*np.pi*e_0*c) # q/(4*np.pi*e_0*c**2)
    A_x1 = schaling_x * A_x1
    A_xn = schaling_x * A_xn
    A_xnXmax = schaling_x * A_xnXmax
    A_xngem = schaling_x * A_xngem




    lines0.append(plt.plot(a,np.abs(A_x1), label= 'dt/dtr')[0])
    lines1.append(plt.plot(a,np.abs(A_xn), label= '1 - n cos(theta)')[0])
    #lines2.append(plt.plot(a,np.abs(A_xngem), label= '1 - <n> cos(theta)')[0])

    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(cti/0.3,a)[0])



plt.legend()
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([0,1e-22])
plt.xlim([20000,42000])
plt.title("Geomagnetic potential against z, s = 1400, theta=85 ")
plt.xlabel("z(m)")
plt.ylabel("Ax")
plt.show()