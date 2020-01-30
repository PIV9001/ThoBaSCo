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
#from crude_sum.py import fullcrude_amp




ec = 1.60217662e-19         #elemental charge
c = 299792458               #speed of light
hbar = 6.582119514e-16      #eV*s
hev = 4.135667662e-15       #Planck Constant h in eV*s
em = 0.51099895e6           #rest mass energy of electreon in eV



#roll bunch here
Ne = 2000#200000       #number of electrons
#gamma0 = 10     #Lorentz factor of mean energy



#filename for reading ASTRA-files
#filename = 'fu1380nm-35-10pc.ini'
#filename = 'fu1380nm-35-10es.ini'
#filename = 'gauss_low_noise.ini'
#filename = 'ful460nm-35-10pc.ini'
filename = 'Sdalinac.ini'
#filename = 'Mesa.ini'


fname = filename.strip('.ini')




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
###             reading ASTRA files                                         ###
###############################################################################


#read ASTRA-file using np.genfromtxt
data = np.genfromtxt(filename)
#print(data[0,:])

#splitting columns in seperate arrays
#coordinates
horpos = data[:, 0].copy()
horpos = np.array(horpos)
#print(np.shape(horpos))
verpos = data[:, 1].copy()
verpos = np.array(verpos)
longpos = data[:, 2].copy()
longpos = np.array(longpos)
#momenta
phor = data[:,3].copy()
phor = np.array(phor)
pvert = data[:,4].copy()
pvert = np.array(pvert)
envec = data[:,5].copy()       #longitudinal momentum
envec = np.array(envec)

# CONVERTING ##################################################################
'''
ASTRA uses the first line as an origin. Coordinates in the first line are
absolute. Coordinates in the following lines are relative to the origin. Here
we convert every line to absolute coordinates.
'''
horpos[0] = 0
verpos[0] = 0
longpos[0] = 0
#
kin = np.copy(envec[0])
#
envec = envec + kin
envec[0] = kin
#convert eV to Lorentz factor
gvec = envec/511e3
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

gamma = gvec[0]
bet0 = np.sqrt(1-1/gamma**2)

###############################################################################
#                       Laser Parameters
###############################################################################
a0 = 2#0.5
nmax = 1        #harmonic
acc = 10
N0 = 7

###   Bill microbunching
#lph = 1380e-9
##l0 = lph*4*gamma**2
#l0 = lph*4*gamma**2
#w0base = 2*np.pi*c/l0
#w0 = a0_freq_corr(w0base, bet0, a0)


####   Mesa
#l0 = 575e-9
##l0 = lph*4*gamma**2
##l0 = lph*4*gvec[0]**2
#w0 = 2*np.pi*c/l0


###   sDalinac
elaser = 1 #photon energy of laser in eV
l0 = hev*c/elaser
#l0 = lph*4*gamma**2
#l0 = lph*4*gvec[0]**2
w0 = 2*np.pi*c/l0


###############################################################################

wmid = 4*gamma**2*w0/(1+a0**2/2)


#gamma = 10
#dgamma = gamma*1e-3
#xvar = 5e-5
#zvar = 1.5e-3
##roll bunch parameters:
##gvec = np.random.normal(gamma, dgamma, Ne)
#gvec = np.ones(Ne)*gamma
#xv = np.random.normal(0, xvar, Ne)
#yv = np.random.normal(0, xvar, Ne)
#zv = np.random.normal(0, zvar, Ne)
##iavec = np.random.normal(incang, dincang, Ne)






###############################################################################
#############   calculation loop  #############################################
###############################################################################

Atsum = np.zeros((xsize,xsize),dtype = np.complex128)
Apsum = np.zeros((xsize,xsize),dtype = np.complex128)
Intsum = np.zeros((xsize,xsize))

dist = 12
#npix = 150


tic = time.perf_counter()

for cnt in range(Ne):#,Ne+1):
    evec = [xv[cnt], yv[cnt], zv[cnt],phor[cnt],pvert[cnt],gvec[cnt]]
    #test case for identical particles
    #evec = [xv[0], yv[0], zv[0],phor[0],pvert[0],gvec[0]]
    tmesh, phi, plmesh, incang2 = pathphase(dist, drange, xsize, evec, lvec, gamma)
    print('incang: ', incang2)
    Atheta, Aphi, w, psi0 = arb_angle_circ_amp(evec, tmesh, phi, a0, w0, gvec[cnt], incang2, N0, nmax)
    #test case for identical particles
    #Atheta, Aphi, w, psi0 = arb_angle_circ_amp(evec, tmesh, phi, a0, w0, gvec[0], incang2, N0, nmax)
    #print(psi0)
#
    #calculate phase based on path length
    plphase = (plmesh/(c/w)%1)*2*np.pi
#    plall[:,:,cnt] = plphase
#    Atsum += Atheta*np.exp(1j*plphase)
#    Apsum += Aphi*np.exp(1j*plphase)
    #incoherent addition
    Intsum += np.absolute(Atheta)**2 + np.absolute(Aphi)**2
    Atsum += Atheta*np.exp(1j*plphase)
    Apsum += Aphi*np.exp(1j*plphase)
    print(cnt)
    
toc = time.perf_counter()
print('elapsed time:')
print(toc-tic)
#
###############################################################################

###############################################################################
#                  1933 plots
###############################################################################

wvec = np.linspace(0.25*wmid,2.5*wmid,150)
tvec = np.linspace(0,1.6,150)/gamma


#fullcrude_amp (w, thetavec, phi, a0, w0, gamma, N0, nmax, acc):
###############################################################################


Int = np.absolute(Atsum)**2 + np.absolute(Apsum)**2
#Int = np.absolute(Apsum + Atsum)**2
#incoherent case:
#Int = Atsum+Apsum

#normalize the intensity:
#normspec = (Int-np.nanmin(Int))/(np.nanmax(Int)-np.nanmin(Int))
#normspec = Int

print('max intensity:')
print(np.max(Int))





###############################################################################
#       saving data to file
###############################################################################
#pre = 'gauss_1ke_'
###np.savetxt(pre+fname+'_Atsum',Atsum)
###np.savetxt(pre+fname+'_Apsum',Apsum)
##np.savetxt(pre+fname+'_coh_Int.txt',Int)
##np.savetxt(pre+fname+'_incoh_Int.txt',Intsum)
#
##np.savetxt(pre+'_coh_Int.txt',Int)
#np.savetxt(pre+'_incoh_Int.txt',Intsum)

#np.savetxt('circ_gauss_1e_incoh_Int.txt',Intsum)


#normspec = np.arctan2(ymesh,xmesh)


###############################################################################
################   plotting    ################################################
###############################################################################
nbin = np.int_(np.sqrt(Ne))   #number o bins in histogram

xysize = 12
zsize = 12


fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax.plot_wireframe(xmesh,ymesh,Int,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'intensity in arb. units',fontsize=zsize)
plt.title('coherent addition',fontsize=zsize)
#
ax = fig.add_subplot(122, projection='3d')
ax.plot_wireframe(xmesh,ymesh,Intsum,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'intensity in arb. units',fontsize=zsize)
plt.title('incoherent addition',fontsize=zsize)





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






