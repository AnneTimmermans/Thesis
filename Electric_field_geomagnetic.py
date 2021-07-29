# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:48:45 2020

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



hoek = range(80,90,2) #[62.98, 71.61, 79.49, 84.95, 87.08]
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
    
    #s = 1000
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
   
    #bepalen dAx/dx 
    avt = (v_0 - v_z)*(1/1225)**alpha * ((alpha*rho_0*C*(R*np.cos(theta) - ctpi))/(np.sqrt(((-ctpi)**2 + R**2 - 2*R*ctpi* np.cos(theta)))) * np.exp(-alpha*C*a)) 
    
    #aNX = Ne*Fp*((np.exp((-(3*X*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0))))/2+X-X_max)/X_0)*(-(3*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0))))/2-(X_0*(X/X_0+(2*X_max)/X_0)*(3/(X_0*(X/X_0+(2*X_max)/X_0))-(3*X)/(X_0**2*(X/X_0+(2*X_max)/X_0)**2)))/2+1))/X_0)
    #aNX = Ne*Fp*(-(np.exp((-(3*X*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0))))/2+X-X_max)/X_0)*((3*X+6*X_max)*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0)))-2*X+2*X_max))/(2*X_0*(X+2*X_max)))
    aNX = Ne*Fp*((np.exp((-(3*X*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0))))/2+X-X_max)/X_0)*((3*X+6*X_max)*np.log((3*X)/(X_0*(X/X_0+(2*X_max)/X_0)))-2*X+2*X_max))/(2*X_0*(X+2*X_max)))
    aXt = rho_0 * np.exp(C*ctpi*np.cos(theta)) * 10**(-4)
    aNt = aNX*aXt
    
    ant = np.exp(C*ctpi*np.cos(theta)) * ((2*k)/(ctpi**3*C*np.cos(theta)) - (2*k)/(ctpi**2) + (k*C*np.cos(theta))/(ctpi)) - (2*k)/(ctpi**3 * C*np.cos(theta))
    asqrtt = ((2*np.cos(theta)**3*(np.cos(theta)*l-z))/np.sin(theta)**2+2*l)/(2*np.sqrt((np.cos(theta)**2*(np.cos(theta)*l-z)**2)/np.sin(theta)**2+l**2))
    aD1t = (ant*np.sqrt(ctpi**2 + s**2) + asqrtt*an)
    #aD1t = -(ant*np.sqrt(ctpi**2 + s**2) + asqrtt*an)
    
    a1_sqrtt = (l*((2*np.cos(theta)**3*(np.cos(theta)*l-z))/np.sin(theta)**2+2*l))/(2*((np.cos(theta)**2*(np.cos(theta)*l-z)**2)/np.sin(theta)**2+l**2)**(3/2))-1/np.sqrt((np.cos(theta)**2*(np.cos(theta)*l-z)**2)/np.sin(theta)**2+l**2)
    aD2t = (an * (ctpi)/(np.sqrt((-ctpi)**2 + s**2)) + a1_sqrtt * n)
    #aD2t = -(an * (ctpi)/(np.sqrt((-ctpi)**2 + s**2)) + a1_sqrtt * n)
    
    aDt = aD1t + aD2t
    
    avN = avt*N + aNt*v
    #avN = 0.4081632653e-3 * (v_z - v_0) * (1 - rho_0 * np.exp(-C * (-R + (ctpi ** 2 + R ** 2 - 2 * R * ctpi * np.cos(theta)) ** 0.5e0)) / 1225) ** alpha * alpha * rho_0 * C * (ctpi ** 2 + R ** 2 - 2 * R * ctpi * np.cos(theta)) ** (-0.5e0) * (2 * ctpi - 2 * R * np.cos(theta)) * np.exp(-C * (-R + (ctpi ** 2 + R ** 2 - 2 * R * ctpi * np.cos(theta)) ** 0.5e0)) / (1 - rho_0 * np.exp(-C * (-R + (ctpi ** 2 + R ** 2 - 2 * R * ctpi * np.cos(theta)) ** 0.5e0)) / 1225) * Ne * Fp * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 + 2 * X_max / X_0))) / X_0) + (v_0 + (v_z - v_0) * (1 - rho_0 * np.exp(-C * (-R + (ctpi ** 2 + R ** 2 - 2 * R * ctpi * np.cos(theta)) ** 0.5e0)) / 1225) ** alpha) * Ne * Fp * (0.73e2 / 0.625e3 * np.exp(C * ctpi * np.cos(theta)) - 0.1752000000e0 * np.exp(C * ctpi * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 + 2 * X_max / X_0)) - 0.5000000000e0 * (0.219e3 / 0.625e3 * np.exp(C * ctpi * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 + 2 * X_max / X_0) - 0.15987e5 / 0.390625e6 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) ** 2 / X_0 ** 2 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 + 2 * X_max / X_0) ** 2) * X_0 * (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 + 2 * X_max / X_0)) / X_0 * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) / X_0 + 2 * X_max / X_0))) / X_0)
    a1_D = (-1/(D**2))*aDt
    #a1_D = -0.1e1 / (k / ctpi ** 2 / C / np.cos(theta) * (1 - np.exp(C * ctpi * np.cos(theta))) + k / ctpi * np.exp(C * ctpi * np.cos(theta)) * (ctpi ** 2 + s ** 2) ** 0.5e0 + (1 - 0.2260000000e-6 * rho_0 / ctpi / C / np.cos(theta) + 0.2260000000e-6 * rho_0 * np.exp(C * ctpi * np.cos(theta)) / ctpi / C / np.cos(theta)) * ctpi * (ctpi ** 2 + s ** 2) ** (-0.5e0) + 1) ** 2 * (-2 * k / ctpi ** 3 / C / np.cos(theta) * (1 - np.exp(C * ctpi * np.cos(theta))) - k / ctpi ** 2 * np.exp(C * ctpi * np.cos(theta)) - k / ctpi ** 2 * np.exp(C * ctpi * np.cos(theta)) * (ctpi ** 2 + s ** 2) ** 0.5e0 + k / ctpi * C * np.cos(theta) * np.exp(C * ctpi * np.cos(theta)) * (ctpi ** 2 + s ** 2) ** 0.5e0 + 0.10e1 * k * np.exp(C * ctpi * np.cos(theta)) * (ctpi ** 2 + s ** 2) ** (-0.5e0) + (0.2260000000e-6 * rho_0 / ctpi ** 2 / C / np.cos(theta) + 0.2260000000e-6 * rho_0 * np.exp(C * ctpi * np.cos(theta)) / ctpi - 0.2260000000e-6 * rho_0 * np.exp(C * ctpi * np.cos(theta)) / ctpi ** 2 / C / np.cos(theta)) * ctpi * (ctpi ** 2 + s ** 2) ** (-0.5e0) + (1 - 0.2260000000e-6 * rho_0 / ctpi / C / np.cos(theta) + 0.2260000000e-6 * rho_0 * np.exp(C * ctpi * np.cos(theta)) / ctpi / C / np.cos(theta)) * (ctpi ** 2 + s ** 2) ** (-0.5e0) - 0.10e1 * (1 - 0.2260000000e-6 * rho_0 / ctpi / C / np.cos(theta) + 0.2260000000e-6 * rho_0 * np.exp(C * ctpi * np.cos(theta)) / ctpi / C / np.cos(theta)) * ctpi ** 2 * (ctpi ** 2 + s ** 2) ** (-0.15e1))
    aAx = avN * (1/D) + a1_D * (v*N)

    
    E_x = aAx
    #E_x = 0.4081632653e-3 * (v_z - v_0) * (1 - rho_0 * np.exp(-C * (-R + (l ** 2 + R ** 2 + 2 * R * l * np.cos(theta)) ** 0.5e0)) / 1225) ** alpha * alpha * rho_0 * C * (l ** 2 + R ** 2 + 2 * R * l * np.cos(theta)) ** (-0.5e0) * (2 * l + 2 * R * np.cos(theta)) * np.exp(-C * (-R + (l ** 2 + R ** 2 + 2 * R * l * np.cos(theta)) ** 0.5e0)) / (1 - rho_0 * np.exp(-C * (-R + (l ** 2 + R ** 2 + 2 * R * l * np.cos(theta)) ** 0.5e0)) / 1225) * Ne * Fp * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0))) / X_0) / ((k / l ** 2 / C / np.cos(theta) * (1 - np.exp(-C * l * np.cos(theta))) - k / l * np.exp(-C * l * np.cos(theta))) * (l ** 2 + s ** 2) ** 0.5e0 - (1 + 0.2260000000e-6 * rho_0 / l / C / np.cos(theta) - 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l / C / np.cos(theta)) * l * (l ** 2 + s ** 2) ** (-0.5e0) + 1) + (v_0 + (v_z - v_0) * (1 - rho_0 * np.exp(-C * (-R + (l ** 2 + R ** 2 + 2 * R * l * np.cos(theta)) ** 0.5e0)) / 1225) ** alpha) * Ne * Fp * (-0.73e2 / 0.625e3 * np.exp(-C * l * np.cos(theta)) + 0.1752000000e0 * np.exp(-C * l * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0)) - 0.5000000000e0 * (-0.219e3 / 0.625e3 * np.exp(-C * l * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0) + 0.15987e5 / 0.390625e6 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) ** 2 / X_0 ** 2 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0) ** 2) * X_0 * (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0)) / X_0 * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0))) / X_0) / ((k / l ** 2 / C / np.cos(theta) * (1 - np.exp(-C * l * np.cos(theta))) - k / l * np.exp(-C * l * np.cos(theta))) * (l ** 2 + s ** 2) ** 0.5e0 - (1 + 0.2260000000e-6 * rho_0 / l / C / np.cos(theta) - 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l / C / np.cos(theta)) * l * (l ** 2 + s ** 2) ** (-0.5e0) + 1) - (v_0 + (v_z - v_0) * (1 - rho_0 * np.exp(-C * (-R + (l ** 2 + R ** 2 + 2 * R * l * np.cos(theta)) ** 0.5e0)) / 1225) ** alpha) * Ne * Fp * np.exp((0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) - X_max - 0.1752000000e0 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) * np.log(0.219e3 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 / (0.73e2 / 0.625e3 / C / np.cos(theta) * np.exp(-C * l * np.cos(theta)) / X_0 + 2 * X_max / X_0))) / X_0) / ((k / l ** 2 / C / np.cos(theta) * (1 - np.exp(-C * l * np.cos(theta))) - k / l * np.exp(-C * l * np.cos(theta))) * (l ** 2 + s ** 2) ** 0.5e0 - (1 + 0.2260000000e-6 * rho_0 / l / C / np.cos(theta) - 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l / C / np.cos(theta)) * l * (l ** 2 + s ** 2) ** (-0.5e0) + 1) ** 2 * ((-2 * k / l ** 3 / C / np.cos(theta) * (1 - np.exp(-C * l * np.cos(theta))) + 2 * k / l ** 2 * np.exp(-C * l * np.cos(theta)) + k / l * C * np.cos(theta) * np.exp(-C * l * np.cos(theta))) * (l ** 2 + s ** 2) ** 0.5e0 + 0.10e1 * (k / l ** 2 / C / np.cos(theta) * (1 - np.exp(-C * l * np.cos(theta))) - k / l * np.exp(-C * l * np.cos(theta))) * (l ** 2 + s ** 2) ** (-0.5e0) * l - (-0.2260000000e-6 * rho_0 / l ** 2 / C / np.cos(theta) + 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l + 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l ** 2 / C / np.cos(theta)) * l * (l ** 2 + s ** 2) ** (-0.5e0) - (1 + 0.2260000000e-6 * rho_0 / l / C / np.cos(theta) - 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l / C / np.cos(theta)) * (l ** 2 + s ** 2) ** (-0.5e0) + 0.10e1 * (1 + 0.2260000000e-6 * rho_0 / l / C / np.cos(theta) - 0.2260000000e-6 * rho_0 * np.exp(-C * l * np.cos(theta)) / l / C / np.cos(theta)) * l ** 2 * (l ** 2 + s ** 2) ** (-0.15e1))
    
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    #e_0 = 8.8541878128 * 10**(-12) # F⋅m−1
    q = 1.6 * 10**(-19)
    #c = 8 * 10**(8)
    schaling = q/(4*np.pi*e_0)
    E_x = schaling * E_x


    
    lines0.append(plt.plot(a,E_x)[0])
    #lines0.append(plt.plot(a,aD)[0])
    #lines0.append(plt.plot(a,atx)[0])



plt.legend(lines0, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([-2e-14,3e-14])
plt.xlim([0,50000])
#plt.xlim([20,140])
#plt.title("E_x against time")
plt.title("E_x against height analytical, zenith = 80")
plt.xlabel("z(m)")
plt.ylabel("E_x")
plt.show()


plt.figure(dpi=100)#Gepasseerde massa
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 80 * (np.pi/180)
    
    delta = 1
    
    h = 1 # m --> hoogte pannenkoek
    
    #s = 1300
    d = s/(np.cos(theta))
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
    ctpi = -l
    ctpi_d = -l + delta
    l_d = -ctpi_d
    d_d = (l_d*np.cos(theta) - z)/(-np.sin(theta))
    s_d = d_d * np.cos(theta)
    
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
    v_d = v_0 + (v_z - v_0)*(1 - (rho_d/1225))**alpha 
    A_x = (v*N)/D
    A_x_d = (v_d * N_d)/D_d
    
    #verhouding Charge excess en geomagnetic
    CG = Cx/v#A_0/A_x
   
    #bepalen dA0/dx 
    E_x = (A_x - A_x_d)/(delta)
    
    aN = (N - N_d)/(delta)
    aD = (D - D_d)/(delta)
    av = (v - v_d)/(delta)
    avN = (v*N - v_d * N_d)/(delta)
    a1_Dt = (1/D - 1/D_d)/delta
    aD1 = (acti_1 - acti_1_d)/delta
    aD2 = (acti_2 - acti_2_d)/delta
    
    A_1 = avN * (1/D) + a1_Dt * (v*N)
    
    e_0 = 55.26349406 # e2⋅GeV−1⋅fm−1
    #e_0 = 8.8541878128 * 10**(-12) # F⋅m−1
    q = 1.6 * 10**(-19)
    #c = 8 * 10**(8)
    schaling = q/(4*np.pi*e_0)
    E_x = schaling * E_x
    


    lines0.append(plt.plot(cti/0.3, E_x)[0])
    #lines0.append(plt.plot(a,aD)[0])
    #lines0.append(plt.plot(a,atx)[0])



plt.legend(lines0, afstanden)
#plt.plot(x,l)
#plt.yscale("log")
#plt.xscale("log")
plt.ylim([-2e-14,3e-14])
#plt.xlim([0,50000])
plt.xlim([20,140])
plt.title("E_x against height numerical, zenith = 80")
plt.xlabel("z(m)")
plt.ylabel("E_x")
plt.show()

