# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:57:12 2023

@author: o.besson
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg
from scipy import stats

#++++++++++ SECTIONS 4 AND 5.1 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
#----- SCENARIO -----
#Number of elements in the array
N = 10
#Inter-element spacing (relative to wavelength)
d = 0.5
#Positions of the antennas
pos = d*np.arange(0,N,1)
pos.shape = (N,1)
#Mainlobe width
theta_3dB = 0.9/(N*d);

#White noise
sigma2 = 1
#Interference
thetaj = np.array([-20,15])/180*np.pi
thetaj.shape = (thetaj.size,1)
INR = np.array([20,20]) #Interference to noise ratio (dB)
Pj = sigma2*(10**(INR/10))
Aj = np.exp(1j*2*np.pi*pos@np.sin(thetaj.T))
#Interference+noise covariance matrix
C = Aj@np.diag(Pj)@np.conj(Aj.T) + sigma2*np.identity(N)
#Signal of interest
theta_s = 0/180*np.pi
SNR = 0 #Signal to noise ratio (dB)
Ps = sigma2*(10**(SNR/10))
a_s = np.exp(1j*2*np.pi*np.sin(theta_s)*pos)
#Total covariance matrix
R = C + Ps*(a_s@np.conj(a_s.T))
#------------------------------------------------------------------------------

#----- AUXILIARY CHANNELS -----------------------------------------------------
#Blocking matrix
a_0 = a_s
B = linalg.null_space(np.conj(a_0.T)) #Blocking matrix N|(N-1) orthogonal to a_0
#Array of DoA values where beampattern is evaluated
tab_theta = np.arange(-90,91,1)/180*np.pi
tab_theta.shape = (tab_theta.size,1)
#Steering matrix (each column is a vector a(theta))
A = np.exp(1j*2*np.pi*pos@np.sin(tab_theta.T))
#Beampatterns of the columns of B
G_B = 20*np.log10(np.abs(np.conj(B.T)@A)).T
fig, ax = plt.subplots()
ax.plot(tab_theta*180/np.pi,G_B,linewidth=2)
ax.set_title(r'Beampatterns of the columns of $\mathbf{B}$')
ax.set_xlabel(r'Angle of arrival $\theta$ (degrees)')
ax.set_ylabel('dB')
ax.set_ylim(-30,10)
for i in range(len(thetaj)):
    ax.axvline(thetaj[i]*180/np.pi,ls='--',color='k')
#Test another matrix B
#Orthogonal (N-1)|(N-1) matrix
U = linalg.orth(stats.norm.rvs(size=[N-1,N-1]) + 1j*stats.norm.rvs(size=[N-1,N-1])) 
B2 = B@U
G_B2 = 20*np.log10(np.abs(np.conj(B2.T)@A)).T
fig, ax = plt.subplots()
ax.plot(tab_theta*180/np.pi,G_B2,linewidth=2)
ax.set_title(r'Beampatterns of the columns of $\mathbf{B}_{2}$')
ax.set_xlabel(r'Angle of arrival $\theta$ (degrees)')
ax.set_ylabel('dB')
ax.set_ylim(-30,10)
for i in range(len(thetaj)):
    ax.axvline(thetaj[i]*180/np.pi,ls='--',color='k')   

#----- CHECK THAT OPTIMAL BEAMFORMER IS PARTIALLY ADAPTIVE --------------------
#MVDR beamformer (known C, a_O=a_s)
a_0 = a_s
w_MVDR = linalg.solve(C,a_0)
w_MVDR = w_MVDR/(np.conj(a_0.T)@w_MVDR)
#Blocking matrix
B = linalg.null_space(np.conj(a_0.T)) #Blocking matrix N|(N-1) orthogonal to a_0
#Check that w_MVDR belongs to subspace spanned by [a_0 , B B^H A_j]
#Orthogonal basis of [a_0 , B B^H A_j] Q^H Q = I_{J+1}
Q = linalg.qr(np.concatenate((a_0,B@np.conj(B.T)@Aj), axis=1),mode='economic')[0] 
print(linalg.norm(np.conj(Q.T)@Q-np.eye(len(thetaj)+1)))
print(linalg.norm(w_MVDR-Q@np.conj(Q.T)@w_MVDR))
#------------------------------------------------------------------------------