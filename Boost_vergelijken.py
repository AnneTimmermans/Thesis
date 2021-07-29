# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:36:12 2021

@author: annee
"""
import  math   as math  
import cmath  as  cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches


#s = 100
#d = 100
z = np.linspace(100000,1,1000)



hoek = range(0,90,15)
afstanden = [150]#range(1000, 2000,200)
lines0 = []
lines1 = []
lines2 = []
lines3 = []
lines4 = []


red_patch = mpatches.Patch(color='red', label= '1 - n(z) cos(theta)')
blue_patch = mpatches.Patch(color='blue', label= 'dt/dtr')


SMALL_SIZE = 8
MEDIUM_SIZE = 13
BIGGER_SIZE = 16

#plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
#plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


plt.figure(dpi=100) #signaal afhankelijk van hoogte
for s in afstanden:
    
    #theta = t * (np.pi/180)
    theta = 45 * (np.pi/180)
    
    d = s/(np.cos(theta))
    #s = d * np.cos(theta)
    
    l = (z - d*np.sin(theta))/(np.cos(theta))
    ctpi = -l
    C = 1.168 * 10**(-4) #m-1
    R = 6300 * 10**3 #m
    rho_0 = 1168 #g/m3
    a = -R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5 #m
    rho = rho_0 * np.exp(-C*a) #meter
    k = 0.226*rho_0*10**(-6) #g/m3
    
    X_max = 700 * 10**4 #g/m2
    L_max = (np.log(rho_0 / (C * X_max * np.cos(theta))))/(C * np.cos(theta)) 
    
    s_max = L_max/np.tan(theta)
    d_max = s_max/(np.cos(theta))
    Z_max = L_max * np.cos(theta) + d_max*np.sin(theta)
    
    n = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (l**2 + R**2 + 2*R*l* np.cos(theta))**0.5)) *10**(-6)
    n_xmax = 1 + 0.226 * rho_0 * np.exp(-C*(-R + (L_max**2 + R**2 + 2*R*L_max* np.cos(theta))**0.5)) *10**(-6)
    n_gem = 1 + (0.226*rho_0*10**(-6))/(l*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*l*np.cos(theta)))/(l*C*np.cos(theta))# 10**(-6) is voor de van meter naar centimeter --> n is dimensie loos
    n_max = 1 + (0.226*rho_0*10**(-6))/(L_max*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*L_max*np.cos(theta)))/(L_max*C*np.cos(theta))
    n_0 = 1 + (0.226*rho_0*10**(-6))/(1*C*np.cos(theta)) - (0.226*rho_0*10**(-6)*np.exp(-C*1*np.cos(theta)))/(1*C*np.cos(theta))
    #s = d * np.cos(theta)

    cti = n_gem * ((-ctpi)**2 + s**2)**0.5 + ctpi

    an = (k/((-ctpi)**2 * C * np.cos(theta))) * (1 - np.exp(C * ctpi * np.cos(theta))) + k/ctpi * np.exp(C * ctpi * np.cos(theta))
    acti_1 = an * ((-ctpi)**2 + s**2)**0.5
    acti_2 = n_gem * (ctpi / (((-ctpi)**2 + s**2)**0.5))
    acti_3 = 1
    acti = acti_1 + acti_2 + acti_3
    
    
    boost_n = 1 - n* (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    boost_nXmax = 1 - n_xmax * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    boost_ngem = 1 - n_gem * (-ctpi)/(np.sqrt(ctpi**2 + s**2))
    
    Diff_boost = np.abs((np.abs(acti) - np.abs(boost_n)))
    Ratio_boost = (np.abs(acti)/np.abs(boost_n)) 


    #print(L_max,Z_max, acti_L)
    #lines0.append(plt.plot(a, Ratio_boost, 'green', label = 'Ratio')[0]) # , label= 'dt/dtr'
    lines1.append(plt.plot(a,np.abs(acti),'blue')[0])
    #lines2.append(plt.plot(a,np.abs(boost_ngem), label= '1 - <n> cos(theta)')[0])
    lines3.append(plt.plot(a,np.abs(boost_n), 'red')[0]) #, label= '1 - n(z) cos(theta)'
    #lines2.append(plt.plot(acti,z)[0])
    #lines3.append(plt.plot(x,l)[0])


plt.legend(handles = [red_patch,blue_patch])
#plt.legend()
#plt.plot(cti/0.3,z) #delen door 0.3 voor nanometer
plt.yscale("log")
#plt.xscale("log")
plt.xlim([0,14000])
#plt.ylim([9e-7,2e-4])
plt.xlabel("z(m)")
plt.ylabel("Boost")
plt.title("Approximations boost against z for s = 150 with theta=45")
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
#plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.show()

