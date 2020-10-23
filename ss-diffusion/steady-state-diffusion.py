#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Diffusion controlled current in a steady-state amperometry experiment

For an unknown reason, running directly this script with ./Butler-Volmer.py generates errors while executing it as a python script with python 3 Bulter-Volmer.py works fine


Description
------------
This program simulates the evolution of the current for a diffusion controlled system

The mass transfer ratio as well as the concentration ratio can be varied

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
title = "Diffusion controlled current for a steady state system"

description = """This program simulates the current for a diffusion controlled system. See Girault equation 7.34

$E=E^\circ +\dfrac{RT}{nF}\ln\left(\dfrac{m_R}{m_O}\\right)+\dfrac{RT}{nF}\ln\left(\dfrac{I_{dc}-I}{I-I_{da}}\\right)$

with $I_{da} = nFAm_RC_R$ and $I_{dc}=- nFAm_OC_O$

and $m_i = \dfrac{D_i}{\delta_i}$

The red dot gives the half wave potential

$E_{1/2}=E^\circ + \dfrac{RT}{nF}\ln\left(\dfrac{m_R}{m_O}\\right)$

The blue dot gives the Nernst potential to help differentiate the effect of both parameters on the curve.

Here, 1 electron exchanged and the standard potential is equal to 0.77 V/ESH. The electrode area is supposed to be equal to 1 m$^2$."""

#===========================================================
# --- Initial parameters  ---------------------
#===========================================================

parameters = {
    'ratm' : widgets.FloatSlider(value=1., description='ratio -- $\dfrac{m_R}{m_O}$', min=0.00001, max=10),
    'ratC' : widgets.FloatSlider(value=1, description='$\dfrac{C_R}{C_O}$', min=0.00001, max=10),
    }
#default parameter for the potential of the couple considered : E0, temperature T, number of exchanged electrons n
E0 = 0.77
T = 298.15
n = 1.
F = constants.physical_constants['Faraday constant'][0]
R = constants.R

#===========================================================
# --- Functions to plot-------------------------------------
#===========================================================
def Ehalf(E0, ratm, ratC, V,n,T,F,R):
    """
    Compute the half potential for the system
    """
    #E_1/2
    Ehalf = E0 + R*T/(n*F) * np.log(ratm)
    return Ehalf


def current(E0, ratm, ratC, V,n,T,F,R):
    """
    Current expression
     ratm is the ratio m_R/m_O
     ratC is the concentration ratio c_O/c_R
     V is the potential
     n the number of exchanged electrons
    """
    A = 1. #Area
    #anodic limit current
    Ida = n * F * A * ratm *ratC 
    #cathodic limit current
    Idc = -n * F * A
    Eh = Ehalf(E0, ratm, ratC, V,n,T,F,R)
    expTerm = np.exp(n*F/(R*T) * (V-Eh))
    return 1/(np.abs(Idc)*(1+expTerm))*(Idc+Ida*expTerm) 

     


#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed 
def plot_data(ratm,ratC):
    lines['i'].set_data(V, current(E0, ratm, ratC, V,n,T,F,R))
    Eh = Ehalf(E0, ratm, ratC, V,n,T,F,R)
    lines['Eha'].set_data(Eh, current(E0, ratm, ratC,Eh,n,T,F,R))
    lines['En'].set_data(E0+R*T/(n*F)*np.log(1/ratC), [0.])
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
lines['i'], = ax.plot([], [], lw=2, color='green')
lines['Eha'], = ax.plot([], [], marker='o', markersize = 3,color='red')
lines['En'], = ax.plot([], [], marker='o', markersize = 3,color='blue')

V = np.linspace(E0-1., E0+1, 1001)
ax.set_xlim(V.min(), V.max())
ax.set_ylim(-1.1, 10.1)

ax.set_xlabel('Tension (V/ESH)')
ax.set_ylabel('Current ratio  $I/I_{dc}$')

param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.07, 0.4, 0.15])
#choose_widget = widgets.make_choose_plot(lines, box=[0.015, 0.25, 0.2, 0.15])
reset_button = widgets.make_reset_button(param_widgets)

if __name__=='__main__':
    plt.show()




