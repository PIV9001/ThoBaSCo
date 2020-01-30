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
from arb_angle_lin_amp import arb_angle_lin_amp2
#from arb_angle_lin_amp import arb_angle_lin_amp
#from arb_angle_circ_amp import pathphase
#from ThoBaSCo_kernel import pathphase_lin
from ThoBaSCo_kernel import pathphase_lin_rotmat
#from crude_sum.py import fullcrude_amp
from ThoBaSCo_kernel import a0_freq_corr
from ThoBaSCo_kernel import rot_x
from ThoBaSCo_kernel import rot_y
from ThoBaSCo_kernel import BFF_calc






ec = 1.60217662e-19         #elemental charge
c = 299792458               #speed of light
hbar = 6.582119514e-16      #eV*s
hev = 4.135667662e-15       #Planck Constant h in eV*s
em = 0.51099895e6           #rest mass energy of electreon in eV

#roll bunch here
Ne = 200000       #number of electrons
#gamma0 = 10     #Lorentz factor of mean energy



#filename for reading ASTRA-files
filename = 'fu1380nm-35-10pc.ini'
#filename = 'fu1380nm-35-10es.ini'
#filename = 'gauss_low_noise.ini'
#filename = 'ful460nm-35-10pc.ini'
#filename = 'Sdalinac.ini'

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
#incangle1: (ThetaI) Angle between direction of laser and direction and
#                    design direction of electron beam
incangle1 = np.pi/2
#incangle2: (nu)  
incangle2 = 0
lvec = np.array([np.sin(incangle1)*np.cos(incangle2), np.sin(incangle1)*np.sin(incangle2), -np.cos(incangle1)])
#lvec = np.array([np.sin(incangle), 0, -np.cos(incangle)])
polvec = np.array([0,1,0])
#linplane: vector inside the plane of laser polarization and perpendicular to
#          direction of laser travel (used for incident angle calculations)
#linplane = np.cross(lvec,polvec)


#print(c/wmid)

xsize = 80
drange = 2.5
xvec = np.linspace(-drange,drange,xsize)
#xvec = np.linspace(-1,1,xsize)*.3
xmesh,ymesh = np.meshgrid(xvec,xvec)



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
#
#
pabsvec = np.sqrt(phor**2+pvert**2+envec**2)
#convert eV to Lorentz factor
#gvec = envec/em
gvec = pabsvec/em
#print(gvec[0:10])
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



###############################################################################
################         debug parameters        ##############################
###############################################################################

#gamma = 35.507e6/511e3
bet0 = np.sqrt(1-1/gamma**2)
#
a0 = 2
#
#gam2 = 35.507e6/511e3
##lph = 1380e-9
#lph = 1380e-9
##l0 = lph*4*gamma**2
#l0 = lph*4*gamma**2
#w0base = 2*np.pi*c/l0
#w0 = a0_freq_corr(w0base, bet0, a0)
#N0 = 7

#elaser = 1 #photon energy of laser in eV
#
#l0 = hev*c/elaser
##l0 = lph*4*gamma**2
##l0 = lph*4*gvec[0]**2
##w0base = 2*np.pi*c/l0
##w0 = a0_freq_corr(w0base, bet0, a0)
#w0 = 2*np.pi*c/l0
#N0 = 7




# interaction angles
angx =  0#-np.pi/4               #rotation about the x-axis
angy =  0#np.pi/3                #rotation about the y-axis

rotx = rot_x(angx)
roty = rot_y(angy)

rotmat = np.dot(roty,rotx)
#rotmat = np.dot(rot_x(-np.pi/4),rot_y(np.pi/4))


# manual offsets for tests
postdx = np.std(horpos)
postdy = np.std(verpos)

#xv = np.ones_like(longpos)*postdx*100
#yv = np.ones_like(longpos)*postdy*100





###############################################################################
#                       Laser Parameters
###############################################################################
#a0 = 2#0.5
#
#gam2 = 35.507e6/em
##lph = 1380e-9
#lph = 1380e-9
##l0 = lph*4*gamma**2
#l0 = lph*4*gam2**2
#w0 = 2*np.pi*c/l0
#N0 = 7

lph = 1380e-9
#l0 = lph*4*gamma**2
l0 = lph*4*gamma**2
w0base = 2*np.pi*c/l0
w0 = a0_freq_corr(w0base, bet0, a0)
N0 = 7

nmax = 1        #harmonic
acc = 10
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


t_el = longpos/c

###############################################################################
#############   calculation loop  #############################################
###############################################################################

#Atsum = np.zeros((xsize,xsize),dtype = np.complex128)
#Apsum = np.zeros((xsize,xsize),dtype = np.complex128)
Intbff = np.zeros((xsize,xsize))
Intincoh = np.zeros((xsize,xsize))

dist = 12
#npix = 150


tic = time.perf_counter()


########################################################
#############   calculating BFF      ###################
########################################################

#def BFF_calc (N0, t_el, warr):

evec = [xv[0], yv[0], zv[0],phor[0],pvert[0],gvec[0]]
tmesh, phi, plmesh, incang1, incang2 = pathphase_lin_rotmat(dist, drange, xsize, evec, lvec, polvec, gamma)
evec = np.array([xv[0], yv[0], zv[0]])
bvec = np.array([phor[0],pvert[0],envec[0]])
Atheta, Aphi, w, psi0 = arb_angle_lin_amp2(evec, bvec, tmesh, phi, a0, w0, gvec[0], N0, acc, nmax)


bff = BFF_calc(Ne, t_el, w)



# back to calculation loop #############################

for cnt in range(1, Ne):
    #evec = [0, 0, 0,0,0,gamma] #debug
    evec = [xv[cnt], yv[cnt], zv[cnt],phor[cnt],pvert[cnt],gvec[cnt]]
#    print(evec)
    #evec = [xv[cnt], yv[cnt], zv[cnt],phort,pvertt,gvec[cnt]]# angle test
##    print(evec)
#    evec = [0, 0, 0,phor[cnt],pvert[cnt],gvec[cnt]]
    tmesh, phi, plmesh, incang1, incang2 = pathphase_lin_rotmat(dist, drange, xsize, evec, lvec, polvec, gamma)
    #angle test
#    incang1 = 0
#    incang2 = np.pi/2
    #    
#    print('incang1: ', incang1)
#    print('incang2: ', incang2)
#    betvec = np.array([[np.sin(incang1)*np.cos(incang2)], [np.sin(incang1)*np.sin(incang2)], [np.cos(incang1)]])
#    betvec = betvec/np.linalg.norm(betvec)
    #(evec, tmesh, phi, a0, w0, gamma, bvec, N0, acc, nmax)
    #Atheta, Aphi, w, psi0 = arb_angle_lin_amp(evec, tmesh, phi, a0, w0, gvec[cnt], betvec, N0, acc, nmax)
#    print(evec)
    evec = np.array([xv[cnt], yv[cnt], zv[cnt]])
    bvec = np.array([phor[cnt],pvert[cnt],envec[cnt]])
    bvec = np.dot(rotmat,bvec)
    Atheta, Aphi, w, psi0 = arb_angle_lin_amp2(evec, bvec, tmesh, phi, a0, w0, gvec[cnt], N0, acc, nmax)
    #print(psi0)
    #
    #incoherent addition
    Inttemp = np.absolute(Atheta)**2 + np.absolute(Aphi)**2
    Intincoh += Inttemp
    #bff addition
    Intbff += Inttemp*(1+bff*Ne)

    print(cnt)
    
toc = time.perf_counter()
print('elapsed time:')
print(toc-tic)
#
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
print(np.max(Intincoh))




###############################################################################
#       saving data to file
###############################################################################

#pre = 'linear_1h_a0_05_'
###np.savetxt(pre+fname+'_Atsum',Atsum)
###np.savetxt(pre+fname+'_Apsum',Apsum)
#np.savetxt(pre+fname+'_coh_Int.txt',Int)
#np.savetxt(pre+fname+'_incoh_Int.txt',Intsum)

#np.savetxt('lin_gauss_1e_incoh_Int.txt',Intsum)

#normspec = np.arctan2(ymesh,xmesh)

path = '/home/Data/HU-Box/Uni/Masterarbeit/python/pictures/05Simulation/bff_microbunching/linear/'
#
#
#
np.save(path+'200ke_lin_intcoh.npy', Intincoh)
np.save(path+'200ke_lin_bff.npy', bff)
np.save(path+'200ke_lin_intbff.npy', Intbff)
np.save(path+'200ke_lin_axis1.npy', xmesh)
np.save(path+'200ke_lin_axis2.npy', ymesh)


###############################################################################
################   plotting    ################################################
###############################################################################
nbin = np.int_(np.sqrt(Ne))   #number o bins in histogram

xysize = 12
zsize = 12


fig = plt.figure()
ax = fig.add_subplot(131, projection='3d')
ax.plot_wireframe(xmesh,ymesh,Intincoh,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'intensity in arb. units',fontsize=zsize)
plt.title('incoherent addition',fontsize=zsize)
#
ax = fig.add_subplot(132, projection='3d')
ax.plot_wireframe(xmesh,ymesh,bff,linewidth=0.3,rstride=1, cstride=1)
ax.set_xlabel(r'$\gamma_0 x^\prime / z^\prime$',fontsize=xysize)
ax.set_ylabel(r'$\gamma_0 y^\prime / z^\prime$',fontsize=xysize)
ax.set_zlabel(r'intensity in arb. units',fontsize=zsize)
plt.title('bunching form factor',fontsize=zsize)
#
ax = fig.add_subplot(133, projection='3d')
ax.plot_wireframe(xmesh,ymesh,Intbff,linewidth=0.3,rstride=1, cstride=1)
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






