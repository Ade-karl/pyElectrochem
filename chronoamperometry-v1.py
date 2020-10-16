#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Chronoamperometry
"""

# Importation des librairies
import numpy as np
import scipy.fftpack as fft
import scipy.integrate as integrate
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.constants as constants
from scipy.interpolate import UnivariateSpline
from scipy.misc import derivative

# Definition des fonctions
def nextC(C,t,D_red,deltat,deltax):
    """
    Compute the concentration for the next time interval t+delta t from the concentration given at t,x
    """
    #The concentration is splined to make it easier to calculate the second derivative
    derC = np.gradient(C[:,t-1],deltax)
    lapC = np.gradient(derC,deltax)
    newC = C[:,t-1] + D_red*deltat*lapC
    newC[0] = 0.
    return newC

def intensity(C,x,t,deltax,deltat,n,A,D_red):
    """
    compute the intensity from the concentration profile
    """
    gradCx,gradCt = np.gradient(C)
    
    F = constants.physical_constants['Faraday constant'][0]
    return n*F*A*D_red * gradCx[0,:]
     
# Programme principal
if __name__ == "__main__":
    Ei = 0 #Initial potential
    Ef = 2 #Potential for sweep
    E0 = 0.77 #Standard potential for the couple
    n = 1 #Number of electrons exchanged
    nu = 50.e-3 #sweep rate as V/s
    D_ox = 6.04e-10 #diffusion coefficient of the oxydant in m^2/s
    D_red = 7.19e-10 #diffusion coefficient of the reductor in m^2/s
    C_ox = 0. #initial concentration of the oxydant at the electrode in mol/L
    C_red = 0.05 #initial concentration of the oxydant at the electrode in mol/L
    A = 1e-4 # Area of the electrode in m^2
    l = 1e-3 #length in meter
    tfin = 100 #ending time in seconds
    samplingx = 100
    samplingt = 1000
    ksi = np.sqrt(D_red/D_ox)
    T = 298.15 #Temperature
    F = constants.physical_constants['Faraday constant'][0]
    R = constants.R
    convertMoll = 1000 #to convert mol/L to mol/m^3

    x,deltax = np.linspace(0,l,samplingx,retstep=True)
    t,deltat = np.linspace(0,tfin,samplingt,retstep=True)

    DM = np.minimum(D_ox,D_red)*deltat/(deltax**2 )
    print('DM : {}'.format(DM))
    if DM > 0.5:
        print("the sampling is too scarce, choose more wisely")
    C = np.ones((samplingx,samplingt)) 
    #Initial condition for C(x,0) : C(x,0) = C_red 
    C[:,0] = C_red* convertMoll * np.ones(samplingx)
    #just after step function for potential
    C[:,1] = C_red*convertMoll * np.ones(samplingx)
    C[0,1] = 0
    #Concentration of ox is null at t>0
    #C[0,1:]=0
    ax1 = plt.subplot(2,2,1)
    ax2 = plt.subplot(2,2,2)
    ax3 = plt.subplot(2,1,2)

    for time in range(2,samplingt):
        derC = np.gradient(C[:,time-1],deltax)
        lapC = np.gradient(derC,deltax)
        C[:,time]=nextC(C,time-1,D_red,deltat,deltax)
        if time >= 2:
            #print('C')
            #print(C[:,time])
            #print('laplacien')
            #print(lapC)
            ax1.plot(x,C[:,time])
            ax2.plot(x,lapC)

    #plt.plot(x,C[:,0])
    #plt.plot(x,C[:,1])
    
    i = intensity(C,x,t,deltax,deltat,n,A,D_red)
    ax3.plot(t,i) 
    plt.show()
    pass

