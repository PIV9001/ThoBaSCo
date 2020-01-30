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
#from spec_dist_amp import spec_dist_amp
from arb_angle_lin_amp import arb_angle_lin_spec
from ThoBaSCo_kernel import a0_freq_corr
from ThoBaSCo_kernel import pathphase_lin_spec
from ThoBaSCo_kernel import rot_x
from ThoBaSCo_kernel import rot_y
#from crude_sum import fullcrude_amp




ec = 1.60217662e-19         #elemental charge
c = 299792458               #speed of light
hbar = 6.582119514e-16      #eV*s
hev = 4.135667662e-15       #Planck Constant h in eV*s
em = 0.51099895e6           #rest mass energy of electreon in eV

#roll bunch here
Ne = 1       #number of electrons

#filename for reading ASTRA-files
filename = 'Mesa.ini'
#filename = 'Sdalinac.ini'
#filename = 'gauss_low_noise.ini'

fname = filename.strip('.ini')

###############################################################################
###             reading ASTRA files                                         ###
###############################################################################
# INPUT #######################################################################
#filename = 'fu1380nm-35-10pc.ini'

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
gvec = envec/em
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

#gamma = gvec[0]
#bet0 = np.sqrt(1-1/gamma**2)


###############################################################################
################         debug parameters        ##############################
###############################################################################

gamma = gvec[0]
bet0 = np.sqrt(1-1/gamma**2)

#a0 = 2

#gam2 = 35.507e6/em
##lph = 1380e-9
#lph = 1380e-9
##l0 = lph*4*gamma**2
#l0 = lph*4*gamma**2
#w0base = 2*np.pi*c/l0
#w0 = a0_freq_corr(w0base, bet0, a0)
#N0 = 7



###############################################################################
###    tests with emittance and energy spread              ####################
###############################################################################
Ne2 = 200
ES = 0#50e-2
emit = (0/180)*np.pi


EScalc = np.std(gvec)/gamma

stdm = 50

epsnk0 = 1e-6
enkvar = epsnk0*stdm


envect = np.linspace((1-ES)*kin,(1+ES)*kin,Ne2)
gvect = envect/em

iavec = np.linspace(0,emit,Ne2)

xvt = np.ones_like(envect)*xv[0]
yvt = np.ones_like(envect)*yv[0]
zvt = np.ones_like(envect)*zv[0]



postdx = np.std(horpos)
postdy = np.std(verpos)
postdz = np.std(longpos)

bfx = postdx**2*gamma/epsnk0

sigx = np.sqrt(enkvar*bfx/gamma)


#xvt = np.random.normal(0,sigx,envect.shape)
yvt = np.random.normal(0,sigx,envect.shape)
#zvt = np.random.normal(0,sigx,envect.shape)

#xvt = np.linspace(-postdx/2,postdx/2,Ne2)*stdm
#yvt = np.linspace(-postdx/2,postdx/2,Ne2)*stdm
#zvt = np.linspace(-postdx/2,postdx/2,Ne2)*stdm

#single electron offset
#xvt = xvt + postdx*stdm
#yvt = yvt + postdy*stdm
#zvt = zvt + postdz*stdm



#xvt = np.ones_like(envect)*postdx*500
#yvt = np.ones_like(envect)*postdy*200

phort = np.ones_like(envect)*phor[0]
pvertt = np.ones_like(envect)*pvert[0]



###############################################################################
###############################################################################

###############################################################################
#                       Laser Parameters
###############################################################################
a0 = 2.


#Mesa
#a0 = .5
#
#
l0 = 575e-9
#l0 = lph*4*gamma**2
#l0 = lph*4*gvec[0]**2
w0 = 2*np.pi*c/l0
N0 = 7

##sDalinac
#a0 = 2
#
#elaser = 1 #photon energy of laser in eV
#
#l0 = hev*c/elaser
##l0 = lph*4*gamma**2
##l0 = lph*4*gvec[0]**2
#w0base = 2*np.pi*c/l0
#w0 = a0_freq_corr(w0base, bet0, a0)
#N0 = 7

#microbunching example
#lph = 1380e-9
#l0 = lph*4*gamma**2
#w0base = 2*np.pi*c/l0
#w0 = a0_freq_corr(w0base, bet0, a0)
#N0 = 7

###############################################################################

#wmid = 4*gamma**2*w0/(1+a0**2/2)
wmid = gamma**2*w0*(1+bet0)**2/(1+.5*a0**2)

#phi = 0

nmax = 3        #harmonic
acc = 10

dist = 12
###############################################################################

wsize = 300

# 9th harmonic
#wvec = np.linspace(0,5,wsize)*(w0*4*gamma**2) #a0 = 0.5
#wvec = np.linspace(0,3.2,wsize)*(w0*4*gamma**2) #a0 = 2

# 3rd harmonic
wvec = np.linspace(0.01,1.2,wsize)*(w0*4*gamma**2) #a0 = 2

#tvec = np.linspace(0,2,wsize)/gamma
tgrange = 3   #Theta*gamma range


#establishing canvas for reference particle
evec0 = np.array([xv[0], yv[0], zv[0]])
bvec0 = np.array([phor[0],pvert[0],envec[0]])
#bvec0 = [envec[0],pvert[0],phor[0]]

#evec = [0, 0, 0,0,0,gamma] #debug
tvplot, phi = pathphase_lin_spec(dist, tgrange, wsize, evec0, bvec0, gamma)



# interaction angles
angx =  0#-np.pi/4                #rotation about the x-axis
angy =  0#np.pi/4                #rotation about the y-axis

rotx = rot_x(angx)
roty = rot_y(angy)

rotmat = np.dot(roty,rotx)


wmesh,tmesh = np.meshgrid(wvec,tvplot)






###############################################################################
#############   calculation loop  #############################################
###############################################################################



#npix = 150


#tic = time.perf_counter()

#old loop

#wall = np.zeros((wsize,wsize,Ne))
#Atall = np.zeros((wsize,wsize,Ne),dtype = np.complex128)
#Apall = np.zeros((wsize,wsize,Ne),dtype = np.complex128)
#plall = np.zeros((wsize,wsize,Ne),dtype = np.complex128)
#
#Atsum = np.zeros((wsize,wsize),dtype = np.complex128)
#Apsum = np.zeros((wsize,wsize),dtype = np.complex128)
#Intsum = np.zeros((wsize,wsize))

#for cnt in range(0, Ne):
#    for nh in range(1,nmax+1):
#        print(nh)
##        tmesh, phi, plmesh = pathphase_1993(dist, drange*dist/gamma, wsize, xv[cnt], yv[cnt], zv[cnt], 0, 0, gvec[cnt], gamma)
#        #evec = [xv[cnt], yv[cnt], zv[cnt],0,0,gvec[cnt]]
#        evec = [0, 0, 0,0,0,gamma]
#        #Atheta, Aphi, w, psi0 = arb_angle_circ_amp(evec, tmesh, phi, a0, w0, gvec[cnt], incang, N0, nmax)
#        Atheta, Aphi, w, psi0 = spec_dist_amp(evec, wvec, tvec, phi, a0, w0, gamma, N0, nh, acc)
#        #print(psi0)
#        #tm2 = spec_dist_amp(evec, wvec, tvec, phi, a0, w0, gamma, N0, nh, acc)
#        Atall[:,:,cnt] = Atheta
#        Apall[:,:,cnt] = Aphi
#        wall[:,:,cnt] = w
#        #calculate phase based on path length
##        plphase = (plmesh/(c/w)%1)*2*np.pi
##        plall[:,:,cnt] = plphase
#        Atsum += Atheta#*np.exp(1j*plphase)
#        Apsum += Aphi#*np.exp(1j*plphase)
#        #Intsum += np.absolute((Apsum + Atsum)**2)
#    print(cnt)
  
#loop for integrals

#Ixsum = np.zeros((wsize,wsize),dtype = np.complex128)
#Izsum = np.zeros((wsize,wsize),dtype = np.complex128)
#Intsum = np.zeros((wsize,wsize))
  
#for cnt in range(0, Ne):
#    #pathlength calculations for electron go here
#    for nh in range(1,nmax+1):
#        print(nh)
##        tmesh, phi, plmesh = pathphase_1993(dist, drange*dist/gamma, wsize, xv[cnt], yv[cnt], zv[cnt], 0, 0, gvec[cnt], gamma)
#        evec = [xv[cnt], yv[cnt], zv[cnt],0,0,gvec[cnt]]
#        #evec = [0, 0, 0,0,0,gamma]
#        #Atheta, Aphi, w, psi0 = arb_angle_circ_amp(evec, tmesh, phi, a0, w0, gvec[cnt], incang, N0, nmax)
#        #arb_angle_lin_spec (evec, wvec, tvec, phi, a0, w0, gamma, N0, acc, nmax):
#        Ix, Iz, kpref = spec_dist_amp(evec, wvec, tvec, phi, a0, w0, gamma, N0, nh, acc)
#        #print(psi0)
#        #tm2 = spec_dist_amp(evec, wvec, tvec, phi, a0, w0, gamma, N0, nh, acc)
#        #calculate phase based on path length
##        plphase = (plmesh/(c/w)%1)*2*np.pi
##        plall[:,:,cnt] = plphase
#        Ixsum += kpref*Ix#*np.exp(1j*plphase)
#        Izsum += kpref*Iz#*np.exp(1j*plphase)
#        #Intsum += np.absolute((Apsum + Atsum)**2)
#    print(cnt)
#    
#toc = time.perf_counter()
#print('elapsed time:')
#print(toc-tic)

###############################################################################
#      arb angle lin calc loop                         ########################
###############################################################################
tic = time.perf_counter()

Intsum = np.zeros((wsize,wsize))

# calculation loop for direct file reads

#for cnt in range(Ne):
#    evec = np.array([xv[cnt], yv[cnt], zv[cnt]])
#    bvec = np.array([phor[cnt],pvert[cnt],envec[cnt]])
#    #
#    tvec, phi = pathphase_lin_spec(dist, tgrange, wsize, evec, bvec, gamma)
#    #
##    print('phi')
##    print(phi)
#    #
#    bvec = np.dot(rotmat,bvec)
#    Atheta, Aphi = arb_angle_lin_spec(evec, bvec, wvec, tvec, phi, a0, w0, gvec[cnt], N0, acc, nmax)
#    #
#    #incoherent addition
#    Intsum += np.absolute(Atheta)**2 + np.absolute(Aphi)**2
#    print(cnt)

# calculation loop for altered beam paraameters

for cnt in range(Ne2):
    evec = np.array([xvt[cnt], yvt[cnt], zvt[cnt]])
    bvec = np.array([phort[cnt],pvertt[cnt],envect[cnt]])
    #
    tvec, phi = pathphase_lin_spec(dist, tgrange, wsize, evec, bvec, gamma)
    #
    bvec = np.dot(rot_y(iavec[cnt]),bvec)
    #
    Atheta, Aphi = arb_angle_lin_spec(evec, bvec, wvec, tvec, phi, a0, w0, gvect[cnt], N0, acc, nmax)
    #
    Intsum += np.absolute(Atheta)**2 + np.absolute(Aphi)**2
    print(cnt)
   
toc = time.perf_counter()
print('elapsed time:')
print(toc-tic)

###############################################################################
###############################################################################



###############################################################################
# direct intensity
###############################################################################

#Int = fullcrude_amp(wvec, tvec, phi, a0, w0, gamma, N0, nmax, acc)



###############################################################################
############          saving data           ###################################
###############################################################################
axis1 = wmesh/(4*gamma**2*w0)
axis2 = tmesh*gamma


#path = '/home/Data/HU-Box/Uni/Masterarbeit/python/pictures/05Simulation/sdalin_spec_offset_comp/' #work pc
path = 'D:/HU-Box/HU-Box/Uni/Masterarbeit/python/pictures/05Simulation/emittance_variation/' # personal desk top

#sdalin:
#name1 = 'sdalin_{}e_'.format(Ne2)
##name2 = 'x_spread_x{}_'.format(stdm)
#name2 = 'y_spread_x{}_'.format(stdm)

#MESA:
name1 = 'MESA_lin_{}e_'.format(Ne2)
#name2 = 'x_emit_x{}_'.format(stdm)
name2 = 'y_emit_x{}_'.format(stdm)

#
np.save(path+name1+name2+'int', Intsum)
np.save(path+name1+name2+'axis1', axis1)
np.save(path+name1+name2+'axis2', axis2)


###############################################################################
################   plotting    ################################################
###############################################################################
axis1 = wmesh/(4*gamma**2*w0)
axis2 = tmesh*gamma


xysize = 12
zsize = 12


#photon energy
#fig = plt.figure()
##plt.rcParams.update({'font.size': 15})
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_wireframe(wvec*hbar/1000,tmesh*gamma,Int_amp,linewidth=0.3,rstride=1, cstride=1)
#ax.set_xlabel(r'photon energy in keV',fontsize=xysize)
#ax.set_ylabel(r'scattering angle $\gamma_0\theta$',fontsize=xysize)
#ax.set_zlabel('intensity',fontsize=zsize)
##plt.title('amplitude',fontsize=zsize)
#
#second plot
#ax = fig.add_subplot(122, projection='3d')
#ax.plot_wireframe(wvec*hbar,tmesh*gamma,Int,linewidth=0.3,rstride=1, cstride=1)
#ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
#ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
#ax.set_zlabel('intensity',fontsize=zsize)

#photon normalized frequency
fig = plt.figure()
#plt.rcParams.update({'font.size': 15})
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(axis1,axis2,Intsum,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'normalized photon frequency $\omega/4\gamma^2\omega_0$',fontsize=xysize)
ax.set_ylabel(r'scattering angle $\gamma_0\theta$',fontsize=xysize)
ax.set_zlabel('intensity',fontsize=zsize)












plt.show()






