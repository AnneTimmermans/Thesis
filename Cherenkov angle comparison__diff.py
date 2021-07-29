# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 23:06:51 2021

@author: annee
"""

import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit

SMALL_SIZE = 8
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

#s = 200
#d = 100
z = np.linspace(100000,1,10000)
#y = np.linspace(-100000,100000, 10000)
#x = np.linspace(-100000,100000, 10000)



hoek = [20, 50, 70]#range(70,90,5)
hoek1 = np.linspace(1,89,10000)
afstanden = [900]#range(100, 2000, 400)
afstanden1 = np.linspace(1, 4000, 10000)
height = np.linspace(100000,1,10000)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []
S_0 = []
cher_hoek = []
cher_hoek2 = []



plt.figure(dpi=100) #signaal afhankelijk van hoogte
for t in hoek1:
    theta = t * (np.pi/180)
    rho_0 = 1168
    C = 1.168 * 10**(-4) #m
    
    X_max = 700 * 10**4
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta))
    
    s = afstanden1
    
    
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
    

    s_max = L_max/np.tan(theta)
    d_max = s_max/(np.cos(theta))
    Z_max = L_max * np.cos(theta) + d_max*np.sin(theta)
    
    
    n_xmax = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (L_max**2 + R**2 + 2*R*L_max* np.cos(theta))**0.5)) *10**(-6)
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    n_0 = 1 + (0.226*rho_0*10**(-6))/(1*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*1*np.cos(theta)))/(1*C*np.cos(theta))

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi
    an_gem = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an_gem * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    
    an_max = (k/((L_max)**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + k/(-L_max) * np.exp(C * (-L_max) * np.cos(theta))
    acti_1_max = an_max * ((L_max)**2 + s**2)**0.5
    acti_2_max = n_max * ((-L_max) / (((L_max)**2 + s**2)**0.5))
    acti_3_max = 1
    acti_max = np.abs(acti_1_max + acti_2_max + acti_3_max)
    
    
    #i = np.where(acti_max == acti_max.min())
    i = np.int(np.asarray(np.where(acti_max == acti_max.min())))
    #print(i)
    S_0.append(s[i])
    l_diff = L_max - d[i]*np.sin(theta)
    Theta_s = np.arctan((s[i]/l_diff)) *(180/np.pi)
    cher_hoek.append(Theta_s) 
    #print(theta_s_plus)
    
    #print("l", l[i], np.sqrt(d[i]**2 - s[i]**2), d[i]*np.sin(theta))
#    print("l2", l_diffi, l_diff)
#    print("Theta", Theta_si, Theta_s)
    
    if theta == 80 *(np.pi/180):
        print(Theta_s, L_max)
    
    
    #lines1.append(plt.plot(z,acti)[0])
    #print(cher_hoek)   



#print(cher_hoek)
#plt.legend() #([lines0, lines1, lines2], ['dtr/dt', 'cos-1(1/n_gem)', 'cos(1/(n_Xmax))'])
#plt.legend(lines1, hoek)
#plt.plot(hoek1, cher_hoek, label = "numerical arctan(s/l)") #delen door 0.3 voor nanometer
#plt.yscale("log")

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
    n_0 = 1 + (0.226*rho_0*10**(-6))/(1*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*1*np.cos(theta)))/(1*C*np.cos(theta))

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
    
    an_0 = (k/((1)**2 * C * np.cos(theta))) * (1 - np.exp(C * (-0) * np.cos(theta))) + k/(-1) * np.exp(C * (-0) * np.cos(theta))
    acti_1_0 = an_0 * ((1)**2 + s**2)**0.5
    acti_2_0 = n_0 * ((-1) / (((1)**2 + s**2)**0.5))
    acti_3_0 = 1
    acti_0 = acti_1_0 + acti_2_0 + acti_3_0

    #s = d * np.cos(theta)
    theta_s_n_xmax = np.arccos(1/n_xmax)* (180/np.pi)
    theta_s_xmax = np.arccos(1/n_max)* (180/np.pi)
    theta_s_max = np.arccos(1/n_gem)* (180/np.pi)
    #theta_s_gem = np.arccos((acti_1 + acti_3)/(n_gem))* (180/np.pi)
    #theta_s_emissie = np.arccos((acti_1_max + acti_3_max)/(n_max))* (180/np.pi)
    #theta_s_emissie = np.nan_to_num(theta_s_emissie)
    
    
    a_s = (-n_max)
    b_s = ((k/((L_max)**2 * C * np.cos(theta))) * (1 - np.exp(C * (-L_max) * np.cos(theta))) + 1)
    c_s = -k * np.exp(C * (-L_max) * np.cos(theta))
    theta_s_min = np.arccos((-b_s - np.sqrt(b_s**2 - 4*a_s*c_s))/(2*a_s)) * (180/np.pi)
    theta_s_plus = np.arccos((-b_s - np.sqrt(b_s**2 + 4*a_s*c_s))/(2*a_s)) * (180/np.pi)
    
    a_s_0 = (-n_0)
    b_s_0 = ((k/((1)**2 * C * np.cos(theta))) * (1 - np.exp(C * (-1) * np.cos(theta))) + 1)
    c_s_0 = -k * np.exp(C * (-1) * np.cos(theta))
    #theta_s_plus_0 = np.arccos((-b_s_0 - np.sqrt(b_s_0**2 + 4*a_s_0*c_s_0))/(2*a_s_0)) * (180/np.pi)
    
    Diff_cher_hoek = np.abs(cher_hoek - theta_s_n_xmax)
    
    print(L_max)
    
    delta_s = np.tan(np.abs(cher_hoek - theta_s_n_xmax)* (np.pi/180)) * L_max
    delta_d = (np.tan(np.abs(cher_hoek - theta_s_n_xmax)* (np.pi/180)) * L_max)/(np.cos(80* (np.pi/180)))
    print(delta_d)

    #lines0.append(plt.plot(np.cos(theta_s_plus),acti_max, label= 'boost tegen theta_s dtr/dt (X_max), ABC (+) (nieuw)')[0])
    #lines0.append(plt.plot(hoek1,theta_s_plus, label= 'dtr/dt (X_max), ABC (+)')[0])
    #lines1.append(plt.plot(hoek1,theta_s_min, label= 'dtr/dt (X_max), ABC (-)')[0])
    #lines2.append(plt.plot(hoek1,theta_s_xmax, label = 'arccos(1/(n_Xmax))')[0])
    #lines3.append(plt.plot(hoek1,theta_s_plus_0, label= 'dtr/dt (0), ABC (+) (nieuw)')[0])
    #lines3.append(plt.plot(hoek1,delta_d, label = 'Difference observer')[0])
    lines4.append(plt.plot(hoek1,delta_s, label = 'Difference s')[0])
    #lines4.append(plt.plot(hoek1,Diff_cher_hoek, label = 'Absolute difference')[0])

#    
    #print(Z_max)
    


plt.legend() #([lines0, lines1, lines2], ['dtr/dt', 'cos-1(1/n_gem)', 'cos(1/(n_Xmax))'])
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim([0,1000])
plt.xlim([0,89])
plt.xlabel("Zenith angle (deg)")
#plt.ylabel("Shift in observer distance (m)")
plt.ylabel("Shift in impact parameter s (m)")
#plt.xlabel("cos(Cherenkov angle)")
#plt.ylabel("Boost op X_max")
plt.title("Zenith angle against shift in impact parameter")
#plt.title('Cherenkov angle against boost op X_max')
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
#plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.show()