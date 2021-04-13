# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 13:31:06 2021

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



hoek = range(0,90,11)
hoek1 = np.linspace(89,1,10000)
afstanden = range(100, 2200, 400)
afstanden1 = np.linspace(1, 2500, 10000)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    s= afstanden1
    
    d = afstanden1/(np.cos(theta))
    #s = d * np.cos(theta)

    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    

    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3

    #s = d * np.cos(theta)
    theta_s_max = np.arccos(1/n_gem) * (180/np.pi)
    theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem)) * (180/np.pi)
    




    lines0.append(plt.plot(afstanden1,theta_s_gem)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines0, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,50])
plt.xlim([0,2500])
plt.xlabel("Distance s")
plt.ylabel("Cherenkov angle")
plt.title("Distance s against Cherenkov angle, complete")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    s= afstanden1
    d = afstanden1/(np.cos(theta))
    #s = d * np.cos(theta)

    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    

    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1.0003 #1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3

    #s = d * np.cos(theta)
    theta_s_max = np.arccos(1/n_gem) * (180/np.pi) * (d/d) 
    theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem)) * (180/np.pi)
    




    lines0.append(plt.plot(afstanden1,theta_s_max)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines0, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,50])
plt.xlim([0,2500])
plt.xlabel("Distance s")
plt.ylabel("Cherenkov angle")
plt.title("Distance s against Cherenkov angle, simplified")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    s= afstanden1
    
    d = afstanden1/(np.cos(theta))
    #s = d * np.cos(theta)

    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    

    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3

    #s = d * np.cos(theta)
    theta_s_max = np.arccos(1/n_gem) * (180/np.pi)
    theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem)) * (180/np.pi)
    
    ratio_theta_s = theta_s_gem/theta_s_max
    




    lines1.append(plt.plot(afstanden1,ratio_theta_s)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines1, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,50])
plt.xlim([0,2500])
plt.xlabel("Distance s")
plt.ylabel("Cherenkov angle")
plt.title("Distance s against ratio Cherenkov angle, complete/simplified")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    s= afstanden1
    d = afstanden1/(np.cos(theta))
    #s = d * np.cos(theta)

    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m
    R = 6300 * 10**3 #m
    rho_0 = 1168
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6)
    

    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3

    #s = d * np.cos(theta)
    theta_s_xmax = np.arccos(1/n_max) * (180/np.pi) * (d/d)
    theta_s_max = np.arccos(1/n_gem) * (180/np.pi) * (d/d) 
    theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem)) * (180/np.pi)
    




    lines0.append(plt.plot(afstanden1,theta_s_xmax)[0])
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(lines0, hoek)
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,50])
plt.xlim([0,2500])
plt.xlabel("Distance s")
plt.ylabel("Cherenkov angle")
plt.title("Distance s against Cherenkov angle, theta = arccos(1/n_Xmax)")
plt.show()
