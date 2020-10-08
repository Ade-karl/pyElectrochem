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
import importlib
from scipy.interpolate import UnivariateSpline
from scipy.misc import derivative
import matplotlib.animation as animation

# Definition des fonctions
def nextC(C,t,D_red,deltat,deltax):
    """
    Compute the concentration for the next time interval t+delta t from the concentration given at t,x
    First the laplacian at step t-1 is taken and the partial derivative equation is propagated for time t as
    C(t)=C(t-1)+D*deltat*lapC(t-1)
    """
    lapC = laplacian(C[:,t-1],deltax) 
    newC = C[:,t-1] + D_red*deltat*lapC
    newC[0] = 0.
    return newC

def intensity(C,x,t,deltax,deltat,n,A,D_red):
    """
    compute the intensity from the concentration profile
    """
    gradCx,gradCt = np.gradient(C,deltax)
    
    F = constants.physical_constants['Faraday constant'][0]
    return n*F*A*D_red * gradCx[0,:]#/deltax

def laplacian(f,deltax):
    """
    Computes the laplacian second derivative with central formula except at borders where forward and backward formulas are used
    """
    out = np.zeros_like(f)
    out[1:-1]=(f[2:]-2*f[1:-1]+f[0:-2])/deltax**2
    out[0]=(f[2]-2*f[1]+f[0])/deltax**2
    out[-1]=(f[-3]-2*f[-2]+f[-1])/deltax**2
    return out
     
    
def Cottrell(t,n,A,C_red,D_red,convertMoll):
    """
    Cottrell law for chronoamperometry i = n*F*A*C_red*sqrt(D_red/(Pi*t))
    """
    F = constants.physical_constants['Faraday constant'][0]
    out = np.zeros_like(t)
    out[1:]= n*F*A*C_red*convertMoll*np.sqrt(D_red/(np.pi*t[1:]))
    out[0] = np.inf
    return out 

def lap(C,deltax,t):
    return laplacian(C[:,t],deltax)
def animate(time,C,lapC,intensity,t,x,convertMoll):
    """
    
    """
    line1.set_data(x, C[:,time]/convertMoll)
    line2.set_data(x, lapC[:,time])
    line3.set_data(t[time], intensity[time])
    ax2.set_ylim(np.min([-1,np.min(lapC[:,time])*1.05]),0.)
    
    return line1,line2,line3
     
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
    tfin = 5 #ending time in seconds
    samplingx = 250
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
    #Initializing the array for the concentrations
    C = np.zeros((samplingx,samplingt)) 
    lapC = np.zeros((samplingx,samplingt)) 

    #Initial condition for C(x,0) : C(x,0) = C_red 
    C[:,0] = C_red* convertMoll * np.ones(samplingx)
    lapC[:,0] = laplacian(C[:,0],deltax) 
    #just after step function for potential
    C[:,1] = C_red*convertMoll * np.ones(samplingx)
    C[0,1] = 0
    lapC[:,1] = laplacian(C[:,1],deltax) 
    #Concentration of ox is null at t>0
    #C[0,1:]=0
    fig, axes = plt.subplots(2,2)
    ax1 = plt.subplot(2,2,1)
    ax2 = plt.subplot(2,2,2)
    ax3 = plt.subplot(2,1,2)

    for time in range(2,samplingt):
        #shape = lapC.shape
        #print('lapC : {}'.format(shape))
        C[:,time]=nextC(C,time,D_red,deltat,deltax)
        lapC[:,time] = laplacian(C[:,time],deltax) 
        if time %10 == 0:
            pass
            #lapC = laplacian(C[:,time-1],deltax) 
            #print('C')
            #print(C[:,time]/convertMoll)
            #print('laplacien')
            #print(lapC)
            #ax1.plot(x,C[:,time]/convertMoll)
            #ax2.plot(x,lapC)

    #plt.plot(x,C[:,0])
    #plt.plot(x,C[:,1])

    #lines to animate
    line1, = ax1.plot([], [] )
    line2, = ax2.plot([], [] )
    line3, = ax3.plot([], [] , marker='o')
    lines = [line1,line2,line3] 

    i_th= Cottrell(t,n,A,C_red,D_red,convertMoll)
    i = intensity(C,x,t,deltax,deltat,n,A,D_red)
   
    

    #labels and legend
    #ax1.legend(loc='lower right')
    ax1.set_xlabel('Distance')
    ax1.set_ylabel('Concentration')
    ax1.set_xlim(0.,l)
    ax1.set_ylim(0.,C_red*1.05)
    #ax2.legend(loc='lower right')
    ax2.set_xlabel('Distance')
    ax2.set_ylabel('Laplacien')
    ax2.set_xlim(0.,l)
    ax3.legend(loc='upper right')
    ax3.plot(t,i,label = 'i_solve') 
    ax3.plot(t,i_th,label='i_theo') 
    ax2.set_xlabel('Time')
    ax2.set_ylabel('currrent')
    ratio = i_th/i
    print('max I_th/I : {}'.format(ratio[1:].max()))
    plt.tight_layout()

    ani = animation.FuncAnimation(fig, animate, fargs =(C,lapC,i,t,x,convertMoll) , frames=range(t.size), interval=deltat/1000,save_count=t.size)
    f = r"chrono.mp4" 
    writermp4 = animation.FFMpegWriter(fps=int(1/deltat)) 
    #ani.save(f, writer=writermp4)



    plt.show()
    pass

