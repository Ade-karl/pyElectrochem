#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants
import matplotlib.animation as animation

def nextC(C,t,D,Cini,E,E0,n,T,alpha,k0,deltat,deltax):
    F = constants.physical_constants['Faraday constant'][0]
    R = constants.R
    newC = np.zeros_like(C[:,0,:])
    for i in range(0,D.size):
        lapC = laplacian(C[:,t-1,i],deltax) 
        newC[:,i] = C[:,t-1,i] + D[i]*deltat*lapC
    theta = n*F/(R*T)*(E[t]-E0)
    k_red = k0 *  np.exp(alpha*theta)
    k_ox  = k0 *  np.exp(-(1-alpha)*theta)
    A = np.array([[1+k_red*deltax/D[0],-k_ox*deltax/D[0] ],[-D[0],-D[1]]])
    B = np.array([C[1,t-1,0],-D[0]*C[1,t-1,0]-D[1]*C[1,t-1,1]])
    sol = np.linalg.solve(A, B)
    sol = sol.reshape(len(sol), 1)
    newC[0] = np.transpose(sol)
    return newC

def intensity(C,x,deltax,deltat,n,A,D):
    gradCx = np.gradient(C[:,:,0],deltax,axis=0)
    
    F = constants.physical_constants['Faraday constant'][0]
    return n*F*A*D[0] * gradCx[0,:]

def laplacian(f,deltax):
    out = np.zeros_like(f)
    out[1:-1]=(f[2:]-2*f[1:-1]+f[0:-2])/deltax**2
    out[0]=(f[2]-2*f[1]+f[0])/deltax**2
    out[-1]=(f[-3]-2*f[-2]+f[-1])/deltax**2
    return out

def lap(C,deltax,t,species):
    return laplacian(C[:,t,species],deltax)


def halfPeriod(Ei,Ef,nu):
    halfPeriod = np.abs((float(Ef)-float(Ei))/nu)
    return halfPeriod
         

def potential(Ei,Ef,nu,samplingt,nCycle):
    halfPer = halfPeriod(Ei,Ef,nu)
    t,deltat = np.linspace(0,halfPer,samplingt,endpoint=False,retstep=True)
    if Ei<Ef:
        Eforward = Ei + nu * t 
        Ebackward = Ef - nu * t 
    elif Ei>Ef:
        Eforward = Ei - nu * t 
        Ebackward = Ef + nu * t 
    else:
        print('The starting and ending potential are the same')
    if nCycle % 0.5 != 0:
        print('The number of cycles should be a multiple of 0.5')
    Ecycle = [list(Eforward),list(Ebackward)]
    Efull = []
    for i in range(0,int(nCycle/0.5)):
        Efull.extend(Ecycle[i%2])
    if nCycle % 1 == 0:
        Efull.append(Ei)
    elif nCycle % 0.5 ==0.:
        Efull.append(Ef)
    return np.array( Efull )

def animate(time,C,intensity,E,t,x,convertMoll):
    Cox.set_data(x, C[:,time,1]/convertMoll)
    Cred.set_data(x, C[:,time,0]/convertMoll)
    Et.set_data(t[time],Efull[time] )
    it.set_data(t[time], intensity[time])
    iE.set_data(E[time], intensity[time])
    return [Cox,Cred,Et,it,iE]
     
if __name__ == "__main__":
    alpha = 0.5 
    k0 = 1e-6 
    Ei = 0. 
    Ef = 1.5 
    E0 = 0.77
    n = 1 
    nu = 50.e-3
    D_ox = 6.04e-10
    D_red = 6.04e-10
    C_ox = 0. 
    C_red = 0.05 
    A = 1e-4 
    l = 5e-4 
    nCycle = 1
    samplingx = 100
    samplingt = 10000 
    T = 298.15 
    saveCsv= True 
    saveMovie = False
    saveNpy = False
    movie = "voltammetry-irrev.mp4"

    F = constants.physical_constants['Faraday constant'][0]
    R = constants.R
    convertMoll = 1000 
    halfPer= halfPeriod(Ei,Ef,nu) 
    x,deltax = np.linspace(0,l,samplingx,retstep=True)
    sizet = int(samplingt*2*nCycle)+1
    t,deltat = np.linspace(0,2*nCycle*halfPer,sizet,retstep=True)
    Efull = potential(Ei,Ef,nu,samplingt,nCycle)
    D = np.array([D_red,D_ox])
    Cini = convertMoll*np.array([C_red,C_ox])
    C = np.zeros((samplingx,sizet,2)) 
    C[:,0,0] = C_red* convertMoll * np.ones(samplingx)
    C[:,0,1] = C_ox* convertMoll * np.ones(samplingx)
    DM = np.minimum(D_ox,D_red)*deltat/(deltax**2 )
    print('DM : {}'.format(DM))
    if DM > 0.5:
        print("the sampling is too scarce, choose more wisely : decrease t sampling or raise x sampling")
    frameselect = int(sizet/(50*t[-1]))
    print('display every {} step, length(s) {}, sizet {} '.format(frameselect,t[-1], sizet))
    print('D*T {} , length^2 {}'.format(np.max(D)*t[-1], l**2))
    print('Psi {}'.format(k0/np.sqrt(np.pi * D_red * n * F * nu / (R*T)) ))

    for time in range(1,sizet):
        C[:,time,:]=nextC(C,time,D,Cini,Efull,E0,n,T,alpha,k0,deltat,deltax)

    i = intensity(C,x,deltax,deltat,n,A,D)
   
    fig, axes = plt.subplots(2,2, figsize =(10,6))
    ax1 = plt.subplot(2,2,1)
    ax1.set_xlabel('Distance')
    ax1.set_ylabel('Concentration')
    ax1.set_xlim(0.,l)
    ax1.set_ylim(0.,C_red*1.05)
    ax2 = plt.subplot(2,2,2)
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('intensity (A)')
    ax2.plot(t,i,label = 'i',color='C0') 
    ax2_2 = ax2.twinx() 
    ax2_2.plot(t,Efull,label = 'E',color='C1') 
    ax2.legend(loc='upper left')
    ax2_2.legend(loc='upper right')
    ax3 = plt.subplot(2,2,3)
    ax3.set_xlabel('Voltage (V)')
    ax3.set_ylabel('intensity (A)')
    ax3.plot(Efull,i,label = 'intensity',marker=None) 
    ax2_3 =plt.subplot(2,2,4)  
    ax2_3.set_xlabel('time (m)')
    ax2_3.set_ylabel('c at electrode')
    ax2_3.plot(t,C[0,:,0]/convertMoll,label = 'red',color='C0') 
    ax2_3.plot(t,C[0,:,1]/convertMoll,label = 'ox',color='C1') 
    ax2_3.legend(loc='upper right')
    Cox, = ax1.plot([], [],color='C2' )
    Cred, = ax1.plot([], [],color='C3' )
    Et, = ax2_2.plot([], [] , marker='o',ms=2)
    it, = ax2.plot([], [] , marker='o',ms=2)
    iE, = ax3.plot([], [] , marker='o',ms=3,color='C0')
    plt.tight_layout()
    ani = animation.FuncAnimation(fig, animate, fargs =(C,i,Efull,t,x,convertMoll), blit=True , frames=range(0,sizet+1,frameselect),interval=20)#,save_count=int(sizet/frameselect))
    writermp4 = animation.FFMpegWriter(fps=50) 
    if saveMovie == True:
        ani.save(movie, writer=writermp4)
    filename ="voltammetry-qr-sweep-{}-k0-{}-alpha-{}-E0-{}-C_ox-{}-C_red-{}-Ei-{}-Ef-{}".format(nu,k0,alpha,E0,C_ox,C_red,Ei,Ef)
    if saveCsv == True :
        np.savetxt(filename+'.csv', np.transpose([Efull,i]), delimiter=",")
    if saveNpy == True:        
        with open(filename+'.npy','wb') as fileOutput:
            np.save(fileOutput, [Efull,i,t])
            np.save(fileOutput, x)
            np.save(fileOutput, C  )
        
    plt.show()
    pass
