# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:24:41 2018

@author: kwl
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from arb_angle_circ_amp import arb_angle_circ_amp
#from arb_angle_circ_amp import pathphase
from ThoBaSCo_kernel import pathphase
from ThoBaSCo_kernel import a0_freq_corr
from ThoBaSCo_kernel import BFF_calc
#from crude_sum.py import fullcrude_amp





ec = 1.60217662e-19         #elemental charge
c = 299792458               #speed of light
em = 0.51099895e6           #rest mass energy of electreon in eV 
hbar = 6.582119514e-16      #eV*s
hev = 4.135667662e-15       #Planck Constant h in eV*s

#roll bunch here
Ne = 1#2400000       #number of electrons
#gamma0 = 10     #Lorentz factor of mean energy



filename = '2019-temp-ref_G8_CSR-angle-4_36p0m.txt'
#filename for reading ASTRA-files
#filename = 'fu1380nm-35-10pc.ini'
#filename = 'fu1380nm-35-10es.ini'
#filename = 'gauss_low_noise.ini'
#filename = 'ful460nm-35-10pc.ini'

#fname = filename.strip('.ini')




#

#1D-Bunch-Form-Factor
#check for correct orientation
#tvec -> column vector
#wvec -> row vec

#if isrow(tvec)
#    tvec = tvec';
#end
#if ~isrow(wvec);
#    wvec = wvec';
#end

#wtmat = tvec*wvec;
#
#bffvec = 1/length(tvec)^2*(sum(sin(wtmat)).^2 + sum(cos(wtmat)).^2);
#bffvec = sqrt(bffvec);
#
#end

#

#6D-Bunch-Form-Factor
#BFF = 
#for cnt in range(0, Ne):
#    BFF += np.exp(-1j*X0[cnt])

#interaction parameters


#bunch-parameters:
incangle = 0#np.pi/2 #incident angle for circular polarization
lvec = np.array([np.sin(incangle), 0, -np.cos(incangle)])


#print(c/wmid)

xsize = 80
drange = 2.5
xvec = np.linspace(-drange,drange,xsize)
#xvec = np.linspace(-1,1,xsize)*.3
xmesh,ymesh = np.meshgrid(xvec,xvec)



###############################################################################
###             reading Opal files                                         ###
###############################################################################

data = np.genfromtxt(filename, skip_header=1)


#splitting columns in seperate arrays
#coordinates
horpos = data[:, 0].copy()
verpos = data[:, 2].copy()
longpos = data[:, 4].copy()
#momenta
phor = data[:,1].copy()
pvert = data[:,3].copy()
plong = data[:,5].copy()       #longitudinal momentum

#gvec = plong*10**6/em
gvec = plong


#postdx = np.std(horpos)
#postdy = np.std(verpos)
#postdz = np.std(longpos)

#
t_el = longpos/c
#
#
xv = horpos
yv = verpos
zv = longpos
#pencil beam:
#xv = np.zeros(Ne)
#yv = np.zeros(Ne)
#zv = np.zeros(Ne)
#phor = np.zeros(Ne)
#pvert = np.zeros(Ne)
#gvec = np.ones(Ne)*kin/511e3
###############################################################################

gamma = np.average(gvec)
bet0 = np.sqrt(1-1/gamma**2)

###############################################################################
#                       Laser Parameters
###############################################################################
a0 = 2

#elaser = 1 #photon energy of laser in eV
#
#l0 = hev*c/elaser
#
#w0base = 2*np.pi*c/l0
#w0 = a0_freq_corr(w0base, bet0, a0)
#N0 = 7


nmax = 1        #harmonic





lph = 63.85e-6
 
l0 = lph*4*gamma**2
w0base = 2*np.pi*c/l0
w0 = a0_freq_corr(w0base, bet0, a0)
N0 = 7

lcorr = c/w0





###############################################################################







###############################################################################
#############   calculation loop  #############################################
###############################################################################

#Atsum = np.zeros((xsize,xsize),dtype = np.complex128)
#Apsum = np.zeros((xsize,xsize),dtype = np.complex128)
Intsum = np.zeros((xsize,xsize))
Int_incoh = np.zeros((xsize,xsize))

dist = 12
#npix = 150


tic = time.perf_counter()

evec = np.array([xv[0], yv[0], zv[0],phor[0],pvert[0],gvec[0]])
tmesh, phi, plmesh, incang2 = pathphase(dist, drange, xsize, evec, lvec, gamma)
print('incang: ', incang2)
Atheta, Aphi, wbff, psi0 = arb_angle_circ_amp(evec, tmesh, phi, a0, w0, gvec[0], incang2, N0, nmax)
#print(psi0)
#
bff = BFF_calc(Ne, t_el, wbff)
#
Inttemp = np.absolute(Atheta)**2 + np.absolute(Aphi)**2
Intsum += Inttemp*(1+bff*Ne)
Int_incoh += np.absolute(Atheta)**2 + np.absolute(Aphi)**2



for cnt in range(1,Ne):#range(Ne):
    evec = np.array([xv[cnt], yv[cnt], zv[cnt],phor[cnt],pvert[cnt],gvec[cnt]])
    tmesh, phi, plmesh, incang2 = pathphase(dist, drange, xsize, evec, lvec, gamma)
    print('incang: ', incang2)
    Atheta, Aphi, w, psi0 = arb_angle_circ_amp(evec, tmesh, phi, a0, w0, gvec[cnt], incang2, N0, nmax)
    #print(psi0)
#
    #calculate phase based on path length
#    plphase = (plmesh/(c/w)%1)*2*np.pi
#    plall[:,:,cnt] = plphase
#    Atsum += Atheta*np.exp(1j*plphase)
#    Apsum += Aphi*np.exp(1j*plphase)
    #incoherent addition
#    bff = BFF_calc(Ne, t_el, w)
    #
    Inttemp = np.absolute(Atheta)**2 + np.absolute(Aphi)**2
    Intsum += Inttemp*(1+bff*Ne)
    Int_incoh += np.absolute(Atheta)**2 + np.absolute(Aphi)**2
    
#    Atsum += Atheta*np.exp(1j*plphase)
#    Apsum += Aphi*np.exp(1j*plphase)
    print(cnt)
    
toc = time.perf_counter()
print('elapsed time:')
print(toc-tic)
#
#Int_bff = Intsum*(1+bff)
###############################################################################

###############################################################################
#                  1933 plots
###############################################################################

#wvec = np.linspace(0.25*wmid,2.5*wmid,150)
#tvec = np.linspace(0,1.6,150)/gamma


#fullcrude_amp (w, thetavec, phi, a0, w0, gamma, N0, nmax, acc):
###############################################################################


#Int = np.absolute(Atsum)**2 + np.absolute(Apsum)**2
#Int = np.absolute(Apsum + Atsum)**2
#incoherent case:
#Int = Atsum+Apsum

#normalize the intensity:
#normspec = (Int-np.nanmin(Int))/(np.nanmax(Int)-np.nanmin(Int))
#normspec = Int

print('max intensity:')
print(np.max(Intsum))





###############################################################################
#       saving data to file
###############################################################################
#pre = 'circular_1h_a_0.5_w_corr'
###np.savetxt(pre+fname+'_Atsum',Atsum)
###np.savetxt(pre+fname+'_Apsum',Apsum)
#np.savetxt('fixed_bERLinpro_a0_2_63.85mum_laser_incoh_Int.txt',Int_incoh)
#np.savetxt('fixed_bERLinpro_a0_2_63.85mum_laser_bff_Int.txt',Intsum)
#np.savetxt('fixed_bERLinpro_a0_2_63.85mum_laser_bff.txt',bff)



#normspec = np.arctan2(ymesh,xmesh)


###############################################################################
################   plotting    ################################################
###############################################################################
nbin = np.int_(np.sqrt(Ne))   #number o bins in histogram

xysize = 12
zsize = 12


fig = plt.figure()
ax = fig.add_subplot(131, projection='3d')
ax.plot_wireframe(xmesh,ymesh,Int_incoh,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'intensity in arb. units',fontsize=zsize)
plt.title('incoherent addition',fontsize=zsize)

#plot bff
ax = fig.add_subplot(132, projection='3d')
ax.plot_wireframe(xmesh,ymesh,bff,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'Bunching Form Factor',fontsize=zsize)
plt.title('Bunching Form Factor',fontsize=zsize)
#
ax = fig.add_subplot(133, projection='3d')
ax.plot_wireframe(xmesh,ymesh,Intsum,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'intensity in arb. units',fontsize=zsize)
plt.title('bff addition',fontsize=zsize)





###############    4 Plots   ###################################################
#fig = plt.figure()
#ax = fig.add_subplot(121, projection='3d')
#ax.plot_wireframe(xmesh,ymesh,normspec,linewidth=0.3,rstride=1, cstride=1)
#ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
#ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
#ax.set_zlabel(r'normalized intensity',fontsize=zsize)
#plt.title('spectrum',fontsize=zsize)
##second plot
##ax2 = fig.add_subplot(132)
##ax2.hist(gvec,nbin)
##ax2.set_xlabel(r'lorentz factor',fontsize=xysize)
###ax2.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
##plt.title('gamma distribution',fontsize=zsize)
##third plot
#ax2 = fig.add_subplot(122, projection='3d')
#ax2.scatter(xv[0:Ne],yv[0:Ne],zv[0:Ne])
#ax2.set_xlabel(r'x position in m',fontsize=xysize)
#ax2.set_ylabel(r'y position in m',fontsize=xysize)
#ax2.set_zlabel(r'z position in m',fontsize=xysize)
#plt.title('electron bunch',fontsize=zsize)
###fourth plot
##ax3 = fig.add_subplot(224)
##ax3.hist(iavec,nbin)
##ax3.set_xlabel(r'incident angle',fontsize=xysize)
###ax3.plt.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
##plt.title('angle distribution',fontsize=zsize)


plt.show()






