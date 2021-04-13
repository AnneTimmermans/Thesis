# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 13:47:16 2020

@author: annee
"""


import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#s = 200
#d = 100
z = np.linspace(100000,1,10000)
#y = np.linspace(-100000,100000, 10000)
#x = np.linspace(-100000,100000, 10000)



hoek = range(80,90,2)
afstanden = range(100, 2200, 400)
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
    
    x = d#np.linspace(-d,d,10000)
    yy = d**2 - x**2
    y = np.sqrt(yy)
    xx = d**2 - y**2
    
#    x = 300
#    y = 400
#    d = np.sqrt(x**2 + y**2)
#    s = d*np.cos(theta)
    
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
    CG = Cx/v#A_0/A_x
   
    #bepalen dA0/dx 
    k2 = 0.15 * (1/1225)**epsilon
    alx = -(x*np.sin(theta))/(d*np.cos(theta))
    asx = (x*np.cos(theta))/d
    c = 3 * 10**8
    atx = (x*np.sin(theta))/(np.cos(theta)*d)
    #atx = (-2*n**3 * s + n)/((2 - 2*(n/c)**2) * np.sqrt(-s**2 + n**2 * s**2 - (cti/c)**2)) * asx
    
    arho_l = ((-rho_0**epsilon * epsilon *C*(l+R*np.cos(theta)))/((l**2 + R**2 + 2*R*l* np.cos(theta))**0.5)) * np.exp(-epsilon*C*a) 
    aCx = k2*arho_l*alx
    
    aXl = -rho_0*np.exp(-C*l*np.cos(theta))*10**(-4)
    aXt = rho_0 * np.exp(C*ctpi*np.cos(theta)) * 10**(-4)
    #aNX = Ne*Fp*(-(np.exp((-(3*X*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0))))/2+X-X_max)/X_0)*((3*X+6*X_max)*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0)))-2*X+2*X_max))/(2*X_0*(X+2*X_max)))
    aNX = Ne*Fp*((np.exp((-(3*X*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0))))/2+X-X_max)/X_0)*((3*X+6*X_max)*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0)))-2*X+2*X_max))/(2*X_0*(X+2*X_max)))
    aNx = aNX*aXl*alx

    
    aCxNx = aCx*N + aNx*Cx
    
    aant = np.exp(C*ctpi*np.cos(theta)) * ((2*k)/(ctpi**3*C*np.cos(theta)) - (2*k)/(ctpi**2) + (k*C*np.cos(theta))/(ctpi)) - (2*k)/(ctpi**3 * C*np.cos(theta))
    asqrtx = (1)/(2*((-ctpi)**2 + s**2)**0.5) *(2*ctpi*atx + 2*s*asx)
    #aD_1x = (((-ctpi)**2 + s**2)**0.5 * aant * atx + an*asqrtx)
    aD_1x = -(((-ctpi)**2 + s**2)**0.5 * aant * atx + an*asqrtx)
    asqrt_1x = atx* (1)/(((-ctpi)**2 + s**2)**(1/2)) + ((-1)/(2*((-ctpi)**2 + s**2)**(3/2)) *(2*ctpi*atx + 2*s*asx)) * ctpi 
    #aD_2x = (an*(ctpi / (((-ctpi)**2 + s**2)**0.5)) * atx + asqrt_1x * n)
    aD_2x = -(an*(ctpi / (((-ctpi)**2 + s**2)**0.5)) * atx + asqrt_1x * n)
    aDx = aD_1x + aD_2x
    a1_Dx = (-1/(D**2))*aDx
    
    aA_0 = aCxNx * (1)/(D) + a1_Dx * (Cx * N)
    E_0 = aA_0
    #E_0 = -0.7500000000e-1 * (rho_0 * np.exp(-C * (-R + ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + R ** 2 + 2 * R * (z - d * np.sin(theta))) ** 0.5e0)) / 1225) ** epsilon * epsilon * C * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + R ** 2 + 2 * R * (z - d * np.sin(theta))) ** (-0.5e0) * (-2 * (z - d * np.sin(theta)) / np.cos(theta) ** 2 * np.sin(theta) - 2 * R * np.sin(theta)) * Ne * Fp * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0))) / X_0) / ((k / (z - d * np.sin(theta)) ** 2 * np.cos(theta) / C * (1 - np.exp(-C * (z - d * np.sin(theta)))) - k / (z - d * np.sin(theta)) * np.cos(theta) * np.exp(C * theta)) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** 0.5e0 - (1 + 0.2260000000e-6 * rho_0 / (z - d * np.sin(theta)) / C - 0.2260000000e-6 * rho_0 * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) / C) * (z - d * np.sin(theta)) / np.cos(theta) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.5e0) + 1) + (0.15e0 + 0.15e0 * (rho_0 * np.exp(-C * (-R + ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + R ** 2 + 2 * R * (z - d * np.sin(theta))) ** 0.5e0)) / 1225) ** epsilon) * Ne * Fp * (0.73e2 / 0.625e3 / np.cos(theta) * np.sin(theta) * np.exp(-C * (z - d * np.sin(theta))) - 0.1752000000e0 / np.cos(theta) * np.sin(theta) * np.exp(-C * (z - d * np.sin(theta))) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0)) - 0.5000000000e0 * (0.219e3 / 0.625e3 / np.cos(theta) * np.sin(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0) - 0.15987e5 / 0.390625e6 / C / np.cos(theta) ** 2 * np.exp(-C * (z - d * np.sin(theta))) ** 2 / X_0 ** 2 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0) ** 2 * np.sin(theta)) * X_0 * (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0)) / X_0 * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0))) / X_0) / ((k / (z - d * np.sin(theta)) ** 2 * np.cos(theta) / C * (1 - np.exp(-C * (z - d * np.sin(theta)))) - k / (z - d * np.sin(theta)) * np.cos(theta) * np.exp(C * theta)) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** 0.5e0 - (1 + 0.2260000000e-6 * rho_0 / (z - d * np.sin(theta)) / C - 0.2260000000e-6 * rho_0 * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) / C) * (z - d * np.sin(theta)) / np.cos(theta) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.5e0) + 1) - (0.15e0 + 0.15e0 * (rho_0 * np.exp(-C * (-R + ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + R ** 2 + 2 * R * (z - d * np.sin(theta))) ** 0.5e0)) / 1225) ** epsilon) * Ne * Fp * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * (z - d * np.sin(theta))) / X_0 + 2 * X_max / X_0))) / X_0) / ((k / (z - d * np.sin(theta)) ** 2 * np.cos(theta) / C * (1 - np.exp(-C * (z - d * np.sin(theta)))) - k / (z - d * np.sin(theta)) * np.cos(theta) * np.exp(C * theta)) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** 0.5e0 - (1 + 0.2260000000e-6 * rho_0 / (z - d * np.sin(theta)) / C - 0.2260000000e-6 * rho_0 * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) / C) * (z - d * np.sin(theta)) / np.cos(theta) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.5e0) + 1) ** 2 * ((2 * k / (z - d * np.sin(theta)) ** 3 * np.cos(theta) / C * (1 - np.exp(-C * (z - d * np.sin(theta)))) * np.sin(theta) - k / (z - d * np.sin(theta)) ** 2 * np.cos(theta) * np.sin(theta) * np.exp(-C * (z - d * np.sin(theta))) - k / (z - d * np.sin(theta)) ** 2 * np.cos(theta) * np.exp(C * theta) * np.sin(theta)) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** 0.5e0 + 0.5e0 * (k / (z - d * np.sin(theta)) ** 2 * np.cos(theta) / C * (1 - np.exp(-C * (z - d * np.sin(theta)))) - k / (z - d * np.sin(theta)) * np.cos(theta) * np.exp(C * theta)) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.5e0) * (-2 * (z - d * np.sin(theta)) / np.cos(theta) ** 2 * np.sin(theta) + 2 * d * np.cos(theta) ** 2) - (0.2260000000e-6 * rho_0 / (z - d * np.sin(theta)) ** 2 / C * np.sin(theta) - 0.2260000000e-6 * rho_0 * np.sin(theta) * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) - 0.2260000000e-6 * rho_0 * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) ** 2 / C * np.sin(theta)) * (z - d * np.sin(theta)) / np.cos(theta) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.5e0) + (1 + 0.2260000000e-6 * rho_0 / (z - d * np.sin(theta)) / C - 0.2260000000e-6 * rho_0 * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) / C) * np.sin(theta) / np.cos(theta) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.5e0) + 0.5e0 * (1 + 0.2260000000e-6 * rho_0 / (z - d * np.sin(theta)) / C - 0.2260000000e-6 * rho_0 * np.exp(-C * (z - d * np.sin(theta))) / (z - d * np.sin(theta)) / C) * (z - d * np.sin(theta)) / np.cos(theta) * ((z - d * np.sin(theta)) ** 2 / np.cos(theta) ** 2 + d ** 2 * np.cos(theta) ** 2) ** (-0.15e1) * (-2 * (z - d * np.sin(theta)) / np.cos(theta) ** 2 * np.sin(theta) + 2 * d * np.cos(theta) ** 2))

    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    #e_0 = 8.8541878128 * 10**(-12) # F⋅m−1
    schaling = 1/(4*np.pi*e_0)
    E_0 = schaling * E_0
    
    
    



    
    lines0.append(plt.plot((cti/0.3),E_0)[0])
    #lines0.append(plt.plot(a,aD)[0])
    #lines0.append(plt.plot(a,atx)[0])



plt.legend(lines0, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([-0.3e7,0.7e7])
plt.xlim([20,140])
#plt.xlim([0,50000])
plt.title("E_0 against time")
plt.xlabel("t(ns)")
plt.ylabel("E_0")
plt.show()


plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    delta = 1
    
    h = 1 # m --> hoogte pannenkoek
    
    d = s/(np.cos(theta))
    d_d = d + delta
    s_d = d_d * np.cos(theta)
    #s = d * np.cos(theta)
    
    
    x = d#np.linspace(-d,d,10000)
    x_d = x + delta
    
    yy = d**2 - x**2
    y = np.sqrt(yy)
    xx = d**2 - y**2
    
#    x = 300
#    y = 400
#    d = np.sqrt(x**2 + y**2)
#    s = d*np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    l_d = (z - d_d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    ctpi_d = -l_d
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    a_d = -R + (l_d**2 + R**2 + 2*R*l_d* np.cos(theta))**0.5
    rho = rho_0 * np.exp(-C*a) #meter
    rho_d = rho_0 * np.exp(-C*a_d) #meter
    k = 0.226*rho_0*10**(-6)
    n = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_d = 1 + (0.226*rho_0*10**(-6))/(l_d*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l_d*np.cos(theta)))/(l_d*C*np.cos(theta))

    #s = d * np.cos(theta)

    cti = n * ((-ctpi)**2 + s**2)**0.5 + ctpi
    cti_d = n_d * ((-ctpi_d)**2 + s_d**2)**0.5 + ctpi_d
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    an_d = (k/((-ctpi_d)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi_d * np.cos(theta))) + k/ctpi_d * np.exp(C * ctpi_d * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_1_d = an_d * ((-ctpi_d)**2 + s_d**2)**0.5
    acti_2 = n * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_2_d = n_d * (ctpi_d / (((-ctpi_d)**2 + s_d**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    acti_d = acti_1_d + acti_2_d + acti_3
    D = acti
    D_d = acti_d
    
    X = (1168/(C * np.cos(theta))) * np.exp(-C * l * np.cos(theta)) * 10**(-4) # macht voor per centimeter
    X_d = (1168/(C * np.cos(theta))) * np.exp(-C * l_d * np.cos(theta)) * 10**(-4)


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
    S_d = (3*X_d/X_0)/(X_d/X_0 + 2 * (X_max/X_0))
    Ft = np.exp((X - X_max - 1.5 * X * np.log(S))/X_0)
    Ft_d = np.exp((X_d - X_max - 1.5 * X_d * np.log(S_d))/X_0)
    b = 1
    L = 3.9 #m
    Fp = (h**b) * np.exp(-2*h/L) * (4/L**2)
    N = Ne * Ft * Fp
    N_d = Ne * Ft_d * Fp
    
    # Charge excess potential
    epsilon = 0.6
    Cx = 0.15+0.15*(rho/1225)**epsilon
    Cx_d = 0.15+0.15*(rho_d/1225)**epsilon
    A_0 = (Cx*N)/D
    A_0_d = (Cx_d*N_d)/D_d
    
    # Geomagnetic potential
    alpha = 1
    v_0 = 0.02 #c
    v_z = 0.04
    v = v_0 + (v_z - v_0)*(1 - (rho/1225))**alpha
    A_x = (v*N)/D
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v#A_0/A_x
   
    #bepalen dA0/dx 
    E_0 = (A_0 - A_0_d)/(delta)
    aN = (N - N_d)/(delta)
    aD = (D - D_d)/(delta)
    aCx = (Cx - Cx_d)/(delta)
    aCxN = (Cx*N - Cx_d * N_d)/delta
    a1_Dx = (1/D - 1/D_d)/delta
    aD1 = (acti_1 - acti_1_d)/delta
    aD2 = (acti_2 - acti_2_d)/delta
    
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    #e_0 = 8.8541878128 * 10**(-12) # F⋅m−1
    schaling = 1/(4*np.pi*e_0)
    E_0 = schaling * E_0



    
    lines0.append(plt.plot((cti/0.3),E_0)[0])
    #lines0.append(plt.plot(a,aD)[0])
    #lines0.append(plt.plot(a,atx)[0])



plt.legend(lines0, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([-0.3e7,0.7e7])
plt.xlim([20,140])
#plt.xlim([0,50000])
plt.title("E_0 against time numerical")
plt.xlabel("t(ns)")
plt.ylabel("E_0")
plt.show()

#aCx = -(3*C*epsilon*np.sin(theta)*rho_0**epsilon*x*(np.sin(theta)*d-z-R*np.cos(theta))*np.exp(-C*epsilon*(np.cos(theta)*np.sqrt((z-np.sin(theta)*d)**2/np.cos(theta)**2+(2*R*(z-np.sin(theta)*d))/np.cos(theta)+R**2)-R)))/(20*1225**epsilon*np.cos(theta)*d*np.sqrt((z-np.sin(theta)*d)**2/np.cos(theta)**2+(2*R*(z-np.sin(theta)*d))/np.cos(theta)+R^2))
#aN = -(73*Ne*Fp*np.sin(theta)*x*np.exp((-(219*np.exp(-C*(z-np.sin(theta)*d))*np.log((219*np.exp(-C*(z-np.sin(theta)*d)))/(625*C*X_0*np.cos(theta)*((73*np.exp(-C*(z-np.sin(theta)*d)))/(625*C*X_0*np.cos(theta))+(2*X_max)/X_0))))/(1250*C*np.cos(theta))+(73*np.exp(-C*(z-np.sin(theta)*d)))/(625*C*np.cos(theta))-X_max)/X_0-C*(z-np.sin(theta)*d))*((3750*C*X_max*np.cos(theta)*np.exp(C*(z-np.sin(theta)*d))+219)*np.log((219*np.exp(-C*(z-np.sin(theta)*d)))/(625*C*X_0*np.cos(theta)*((73*np.exp(-C*(z-np.sin(theta)*d)))/(625*C*X_0*np.cos(theta))+(2*X_max)/X_0)))+1250*C*X_max*np.cos(theta)*np.exp(C*(z-np.sin(theta)*d))-146))/(1250*X_0*np.cos(theta)*d*(1250*C*X_max*np.cos(theta)*np.exp(C*(z-np.sin(theta)*d))+73)) 
#aD1 = (2*k*np.cos(theta)*np.sin(theta)*x*(1-np.exp(-C*(z-np.sin(theta)*d))))/(C*d*(z-np.sin(theta)*d)**3)-(k*np.cos(theta)*np.sin(theta)*x*np.exp(-C*(z-np.sin(theta)*d)))/(d*(z-np.sin(theta)*d)**2)
#aD2 = (2*k*np.sin(theta)*x*(1-np.exp(-(k*np.exp(-(C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))**2)/np.sqrt((z-(np.sin(theta)*d)/np.cos(theta))**2+np.cos(theta)**2*d**2)))/(z-(np.sin(theta)*d)/np.cos(theta))-C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta)))))/(C*np.cos(theta)**2*d*(z-(np.sin(theta)*d)/np.cos(theta))**3)-(k*(-(k*((2*C*np.sin(theta)*x*(z-(np.sin(theta)*d)/np.cos(theta)))/(d*np.sqrt((z-(np.sin(theta)*d)/np.cos(theta))**2+np.cos(theta)**2*d**2))+(C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))**2*(2*np.cos(theta)**2*x-(2*np.sin(theta)*x*(z-(np.sin(theta)*d)/np.cos(theta)))/(np.cos(theta)*d)))/(2*((z-(np.sin(theta)*d)/np.cos(theta))**2+np.cos(theta)**2*d**2)**(3/2)))*np.exp(-(C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))**2)/np.sqrt((z-(np.sin(theta)*d)/np.cos(theta))**2+np.cos(theta)**2*d**2)))/(z-(np.sin(theta)*d)/np.cos(theta))-(k*np.sin(theta)*x*np.exp(-(C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))**2)/np.sqrt((z-(np.sin(theta)*d)/np.cos(theta))**2+np.cos(theta)**2*d**2)))/(np.cos(theta)*d*(z-(np.sin(theta)*d)/np.cos(theta))**2)+(C*np.sin(theta)*x)/d)*np.exp(-(k*np.exp(-(C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))**2)/np.sqrt((z-(np.sin(theta)*d)/np.cos(theta))**2+np.cos(theta)**2*d**2)))/(z-(np.sin(theta)*d)/np.cos(theta))-C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))))/(C*np.cos(theta)*(z-(np.sin(theta)*d)/np.cos(theta))**2)
#aD = (-1)/(D**2) * (aD1 + aD2)
#aCxN = aCx*N + aN*Cx
#aA_0 = aCxN * (1)/(D) + aD * (Cx * N)
#E_0 = aA_0
