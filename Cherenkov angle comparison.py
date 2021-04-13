# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:26:55 2021

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



hoek = [80]#range(70,90,5)
hoek1 = np.linspace(89,1,10000)
afstanden = [1300]#range(100, 2000, 400)
afstanden1 = np.linspace(1, 2400, 10000)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for s in afstanden:
    
    theta = hoek1 * (np.pi/180)
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
    

    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    s_max = L_max/np.tan(theta)
    
    n_xmax = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (L_max**2 + R**2 + 2*R*L_max* np.cos(theta))**0.5)) *10**(-6)
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    
    an_max = (k/((L_max)**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_max = an_max * ((L_max)**2 + s_max**2)**0.5
    acti_2_max = n_max * ((-L_max) / (((L_max)**2 + s_max**2)**0.5))
    acti_3_max = 1
    acti_max = acti_1_max + acti_2_max + acti_3_max

    #s = d * np.cos(theta)
    theta_s_n_xmax = np.arccos(1/n_xmax)* (180/np.pi)
    theta_s_xmax = np.arccos(1/n_max)* (180/np.pi)
    theta_s_max = np.arccos(1/n_gem)* (180/np.pi)
    theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem))* (180/np.pi)
    theta_s_emissie = np.arccos((acti_1_max + acti_3_max)/(n_max))* (180/np.pi)
    print(theta_s_gem)
    




    #lines0.append(plt.plot(hoek1,theta_s_gem, label = 'dtr/dt = 0')[0])
    #lines1.append(plt.plot(hoek1,theta_s_max, label= 'arccos(1/n_gem)')[0])
    lines2.append(plt.plot(hoek1,theta_s_xmax, label = 'arccos(1/(n_Xmax))')[0])
    lines3.append(plt.plot(hoek1,theta_s_emissie, label = 'dtr/dt (X_max) = 0')[0])
    lines4.append(plt.plot(hoek1,theta_s_n_xmax, label = 'arccos(1/(n(Xmax)))')[0])
    


plt.legend() #([lines0, lines1, lines2], ['dtr/dt', 'cos-1(1/n_gem)', 'cos(1/(n_Xmax))'])
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,50])
plt.xlim([0,90])
plt.xlabel("Zenith angle")
plt.ylabel("Cherenkov angle")
plt.title("Zenith angle against Cherenkov angle, comparison for different calculations (s = 1300)")
plt.show()

plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek:
    
    theta = t * (np.pi/180)
    #theta = 80 * (np.pi/180)
    
    s = 1300#afstanden
    
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
    

    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    Z_max = L_max * np.cos(theta) + d*np.sin(theta)
    s_max = L_max/np.tan(theta)
    
    n_xmax = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (L_max**2 + R**2 + 2*R*L_max* np.cos(theta))**0.5)) *10**(-6)
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    
    an_max = (k/((L_max)**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_max = an_max * ((L_max)**2 + s_max**2)**0.5
    acti_2_max = n_max * ((-L_max) / (((L_max)**2 + s_max**2)**0.5))
    acti_3_max = 1
    acti_max = acti_1_max + acti_2_max + acti_3_max

    #s = d * np.cos(theta)
    theta_s_n_xmax = np.arccos(1/n_xmax)* (180/np.pi) *(z/z) 
    theta_s_xmax = np.arccos(1/n_max)* (180/np.pi) * (z/z)
    theta_s_max = np.arccos(1/n_gem)* (180/np.pi) 
    theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem))* (180/np.pi) 
    theta_s_emissie = np.arccos((acti_1_max + acti_3_max)/(n_max))* (180/np.pi) *(z/z)
    print(Z_max)
    




    lines0.append(plt.plot(a,theta_s_gem, label = 'dtr/dt = 0')[0])
    lines1.append(plt.plot(a,theta_s_max, label= 'arccos(1/n_gem)')[0])
    #lines2.append(plt.plot(z,theta_s_xmax, label = 'arccos(1/(n_Xmax))')[0])
    #lines3.append(plt.plot(z,theta_s_emissie, label = 'dtr/dt (X_max) = 0')[0])
    #lines4.append(plt.plot(z,theta_s_n_xmax, label = 'arccos(1/(n(Xmax)))')[0])


plt.legend() #([lines0, lines1, lines2], ['dtr/dt', 'cos-1(1/n_gem)', 'cos(1/(n_Xmax))'])
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,50])
plt.xlim([0,50000])
plt.xlabel("Height z")
plt.ylabel("Cherenkov angle")
plt.title("Zenith angle against Cherenkov angle, comparison for different calculations (theta = 80, s = 1300)")
plt.show()
