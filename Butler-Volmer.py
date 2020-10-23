#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Butler-Volmer relationship

For an unknown reason, running directly this script with ./Butler-Volmer.py generates errors while executing it as a python script with python 3 Bulter-Volmer.py works fine


Description
------------
This program simulates the evolution of the Butler-Volmer relationship for a rapid system where the system is always limited by the electron transfer (ie c(0) = C(infty) )

The exchange current as well as the transfer coefficient can be varied and the different curves can be selected.

The convention adopted is the IUPAC one (i>O for an oxydation)


Informations
------------
Author : Martin Vérot  from the ENS de Lyon, France, with the help of some scripts to create the buttons and sliders taken from https://github.com/araoux/python_agregation (written by P Cladé, A Raoux and F Levrier)
Licence : Creative Commons CC-BY-NC-SA 4.0 

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
title = "Butler-Volmer relationship"

description = """This program simulates the i-E curves when the current is limited by the electronic transfer. 

$I=I_0\left( \exp\left( \dfrac{\\alpha n F \eta}{RT} \\right)  - \exp\left( -\dfrac{\left( 1-\\alpha \\right) n F \eta}{RT} \\right)\\right)$

$\eta=E-E^0$, $\\alpha$ is the transfer coefficient for the oxydation

Here, 1 electron exchanged and the standard potential is equal to 0.77 V/ESH. The electrode area is supposed to be equal to 1 m$^2$."""

#===========================================================
# --- Initial parameters  ---------------------
#===========================================================

parameters = {
    'alpha' : widgets.FloatSlider(value=0.5, description='transfer coefficient -- $\\alpha$', min=0, max=1),
    'logI' : widgets.FloatSlider(value=-1.35, description='$\log(I_0)$', min=-5, max=5),
    }
#default parameter for the potential of the couple considered : E0, temperature T, number of exchanged electrons n
E0 = 0.77
T = 298.15
n = 1
F = constants.physical_constants['Faraday constant'][0]
R = constants.R

#===========================================================
# --- Functions to plot-------------------------------------
#===========================================================

def currentOx(E0, alpha, logI, V,n,T,F,R):
    """
    Oxidation current as a function of the voltage, the input for I° is transformed from a logarithmic scale to it corresponding value first
    """
    I = 10**logI
    return  I*np.exp(alpha*n*F*(V-E0)/(R*T)) 

def currentRed(E0, alpha, logI, V,n,T,F,R):
    """
    Reduction current as a function of the voltage, the input for I° is transformed from a logarithmic scale to it corresponding value first
    """
    I = 10**logI
    return  -I*np.exp(-(1-alpha)*n*F*(V-E0)/(R*T)) 

def currentTotal(E0, alpha, logI, V,n,T,F,R):
    """
    Total current (sum of oxydation and reduction current)
    """
    return currentOx(E0, alpha, logI, V,n,T,F,R)+currentRed(E0, alpha, logI, V,n,T,F,R)
     


#===========================================================
# --- Plot of the updated curves ---------------------------
#===========================================================

# This function is called when the sliders are changed 
def plot_data(alpha, logI):
    lines['i-ox'].set_data(V,currentOx(E0, alpha, logI, V,n,T,F,R))
    lines['i-red'].set_data(V,currentRed(E0, alpha, logI, V,n,T,F,R))
    lines['i-tot'].set_data(V,currentTotal(E0, alpha, logI, V,n,T,F,R))
    fig.canvas.draw_idle()


#===========================================================
# --- Initialization of the plot ---------------------------
#===========================================================

fig = plt.figure(figsize=(12,6))
fig.suptitle(title)
#plot of the text
fig.text(0.01, .9, widgets.justify(description), multialignment='left', verticalalignment='top')
#fig.text(0.02, 0.6,"$I=I_0\left( \exp\left( \dfrac{\\alpha n F \eta}{RT} \\right)  - \exp\left( -\dfrac{\left( 1-\\alpha \\right) n F \eta}{RT} \\right)\\right)$")
#fig.text(0.02, 0.6,"$I=I_0\left( \exp\left( \dfrac{\\alpha n F \eta}{RT} \\right)  - \exp\left( -\dfrac{\left( 1-\\alpha \\right) n F \eta}{RT} \\right)\\right)$")

ax = fig.add_axes([0.35, 0.3, 0.6, 0.6])
ax.axhline(0, color='k')

lines = {}
lines['i-ox'], = ax.plot([], [], linestyle='--', lw=2, color='red', visible=True)
lines['i-red'], = ax.plot([], [], linestyle='--', lw=2,color='blue', visible=True)
lines['i-tot'], = ax.plot([], [], lw=2, color='green')

V = np.linspace(E0-1., E0+1, 1001)
ax.set_xlim(V.min(), V.max())
ax.set_ylim(-1.1, 1.1)

ax.set_xlabel('Tension (V/ESH)')
ax.set_ylabel('Current (A)')

param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.07, 0.4, 0.15])
choose_widget = widgets.make_choose_plot(lines, box=[0.015, 0.25, 0.2, 0.15])
reset_button = widgets.make_reset_button(param_widgets)

if __name__=='__main__':
    plt.show()




