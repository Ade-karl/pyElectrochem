#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Diffusion controlled current in a steady-state amperometry experiment

For an unknown reason, running directly this script with ./Butler-Volmer.py generates errors while executing it as a python script with python 3 Bulter-Volmer.py works fine


Description
------------
This program simulates the evolution of the current for a system where both the diffusion and the kinetic regime apply 

The diffusion coefficient ratio,  the concentration ratio, the transfer coefficient and the current exhange can be varied

The convention adopted is the IUPAC one (i>O for an oxydation)


Informations
------------
Author : Martin Vérot  from the ENS de Lyon, France, with the help of some scripts to create the buttons and sliders taken from https://github.com/araoux/python_agregation (written by P Cladé, A Raoux and F Levrier)
Licence : Creative Commons CC-BY-NC-SA 4.0 

Thanks to V Wieczny, R Grüber and J Galiana for some fruitful discussions. The equations used are those from Girault "Analytical and physical electrochemistry" 

WARNING this program requires the widgets.py file to work
"""

import matplotlib.pyplot as plt
import numpy as np
import widgets
import scipy.constants as constants
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
title = "Current for a system under both diffusion and kinetic control"

description = """See equations 7.12, 7.13, 7.59

$\dfrac{1}{I}=\dfrac{1}{nFA\left(k_aC_R(\infty)-k_cC_O(\infty)\\right)}\left( 1 + \dfrac{k_aD_R^{-2/3}+k_cD_O^{-2/3}}{0.62\\nu^{-1/6}\omega^{1/2}}  \\right)$

Here, 1 electron exchanged and the standard potential is equal to 0.36 V/ESH (ferri/ferrocyanide). The electrode area is supposed to be equal to 1 m$^2$ and the system is on a rotating disk electrode in water."""

#===========================================================
# --- Initial parameters  ---------------------
#===========================================================

parameters = {
    'ratD' : widgets.FloatSlider(value=1., description='ratio -- $\dfrac{D_R}{D_O}$', min=0.00001, max=10),
    'ratC' : widgets.FloatSlider(value=1, description='$\dfrac{C_R}{C_O}$', min=0.00001, max=10),
    'alpha' : widgets.FloatSlider(value=0.5, description='transfer coefficient -- $\\alpha$', min=0, max=1),
    'k0' : widgets.FloatSlider(value=-1.35, description='$\log(k^\circ)$', min=-10, max=-1),
    }
#default parameter for the potential of the couple considered : E0, temperature T, number of exchanged electrons n
E0 = 0.36
T = 298.15
n =1
A = 1. #Area
nu = 0.8927e-6 #kinematic viscosity of water
omega = 500.*2*np.pi/60. #rotation speed of 500rpm
Dox = 0.726e-9  #diffusion coeffienct of ferricyanide
C0 = 0.1 #concentration of the oxidant= 1.
F = constants.physical_constants['Faraday constant'][0]
R = constants.R

#===========================================================
# --- Functions to plot-------------------------------------
#===========================================================
def kanodicreduced(E0,V,k0, alpha,n,F,R,T):
    """
    anodic kinetic constant
    """
    return np.exp(alpha * n*F*(V-E0)/(R*T))
def kcathodicreduced(E0,V,k0, alpha,n,F,R,T):
    """
    cathodic kinetic constant
    """
    return np.exp(-(1-alpha) * n*F*(V-E0)/(R*T))
        
def layer(omega, Dox, nu):
    """
    thickness of the diffusion layer where 
        the diffusion coefficient is given in m^2/s , 
        as the kinematic viscosity 
        omega is given in rad/s
    """
    return 1.61 * Dox**(1./3.)*omega**(-1./2.)*nu**(1./6.)

def IcathodicMax(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0, A):
    """
    maximum cathodic current
    """
    fk0 = 10**k0
    delta = layer(omega,Dox,nu)
    Icmax= n * F * A * Dox * C0 / delta
    return Icmax
         
         

def current(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0, A):
    """
    Current expression
     ratm is the ratio m_R/m_O
     ratC is the concentration ratio c_O/c_R
     V is the potential
     n the number of exchanged electrons
    """
    fk0 = 10.**k0
    delta = layer(omega,Dox,nu)
    Icmax=IcathodicMax(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0, A) 
    ka =fk0 * kanodicreduced(E0,V,k0, alpha,n,F,R,T) 
    kc =fk0 * kcathodicreduced(E0,V,k0, alpha,n,F,R,T) 
    kineticTerm =  n*F*A* C0*  (ka*ratC- kc )     
    diffusionTerm = 1. + delta /Dox * (ka /ratD +kc )
    I = kineticTerm /diffusionTerm
    return I/Icmax



def currentButlerVolmer(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0,A):
    """
    Expression of the current for the Butler-VOlmer equation
    """
    #cathodic limit current
    Icmax=IcathodicMax(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0, A) 
    fk0 = 10**k0
    ka = kanodicreduced(E0,V,k0, alpha,n,F,R,T) 
    kc = kcathodicreduced(E0,V,k0, alpha,n,F,R,T) 
    i=n * F * A * fk0 *C0* (ratC * ka - kc ) 
    return i/Icmax

def Ehalf(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0,A):
    """
    Compute the half potential for the system
    """
    #E_1/2
    Ehalf = E0 + R*T/(n*F) * np.log(ratD)
    return Ehalf


def currentDiff(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0,A):
    """
    Current expression
     ratD is the ratio D_R/D_O
     ratC is the concentration ratio c_O/c_R
     V is the potential
     n the number of exchanged electrons
    """
    #cathodic limit current
    Icmax=IcathodicMax(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0, A) 
    
    delta = layer(omega,Dox,nu)
    #anodic limit current
    Ida = n * F * A * Dox*  ratD *ratC*C0/delta 
    #cathodic limit current
    Idc = -n * F * A* Dox* C0/delta
    Eh = Ehalf(E0, ratD, ratC,alpha,k0, V,n,T,F,R,nu,omega,Dox,C0,A)
    expTerm = np.exp(n*F/(R*T) * (V-Eh))
    return 1/(Icmax*(1+expTerm))*(Idc+Ida*expTerm) 


     
     

#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed 
def plot_data(alpha,k0,ratD,ratC):
    lines['i-BV'].set_data(V,currentButlerVolmer(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0,A) )
    lines['i-Diff'].set_data(V,currentDiff(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0,A) )
    lines['i'].set_data(V,current(E0, ratD, ratC, alpha, k0, V,n,T,F,R,nu,omega,Dox,C0,A) )
    fig.canvas.draw_idle()


#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================

fig = plt.figure(figsize=(12,6))
fig.suptitle(title)
#plot of the text
fig.text(0.01, .9, widgets.justify(description), multialignment='left', verticalalignment='top')

ax = fig.add_axes([0.35, 0.3, 0.6, 0.6])
ax.axhline(0, color='k')

lines = {}
lines['i-BV'], = ax.plot([], [], lw=2, color='red')
lines['i-Diff'], = ax.plot([], [], lw=2, color='blue')
lines['i'], = ax.plot([], [], lw=2, color='green')

V = np.linspace(E0-1., E0+1, 1001)
ax.set_xlim(V.min(), V.max())
ax.set_ylim(-10.1, 10.1)

ax.set_xlabel('Tension (V/ESH)')
ax.set_ylabel('Current ratio  $I/|I_{dc}|$')

param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.07, 0.4, 0.15])
choose_widget = widgets.make_choose_plot(lines, box=[0.015, 0.25, 0.2, 0.15])
reset_button = widgets.make_reset_button(param_widgets)

if __name__=='__main__':
    plt.show()




