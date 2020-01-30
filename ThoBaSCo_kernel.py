# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 12:25:21 2018

@author: Paul Volz

ThoBaSCo stands for Thomson BackScattering Code

This Kernel contains all the functions used by this code.

"""
import numpy as np

em = 0.51099895e6           #rest mass energy of electreon in eV

def rot_x(phi):# 3D rotation matrix. rotation about the x-axis
    return np.array([[1.,0.,0.],[0.,np.cos(phi),-np.sin(phi)],[0.,np.sin(phi),np.cos(phi)]])

def rot_y(phi):# 3D rotation matrix. rotation about the y-axis
    return np.array([[np.cos(phi),0.,np.sin(phi)],[0.,1.,0.],[-np.sin(phi),0.,np.cos(phi)]])


def co_trafo (thetaI, tmesh, phi, nu):
    '''
    basic coordinate transformation as in 1995 paper
    
    thetaI - incidence angle, angle between laser and electron beam
             angle between z and z'
    nu -     incidence angle, angle between laser and electron beam
             angle between x and x' (circ pol: nu = thetaI)
    tmesh - polar angle (detector)
    phi - azimuth angle (detector)
    '''
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    #sti = np.sin(thetaI)# not needed
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    snu = np.sin(nu)
    cnu = np.cos(nu)
    #
    l1 = cnu
    l2 = -snu*np.sqrt(1-cti**2/cnu**2)
    l3 = np.tan(nu)*cti
    m1 = 0
    m2 = cti/cnu
    m3 = np.sqrt(1-cti**2/cnu**2)
    n1 = -snu
    n2 = -cnu*np.sqrt(1-cti**2/cnu**2)
    n3 = cti
    #
    Omx = l1*sth*cph+l2*sth*sph+l3*cth
    Omy = m1*sth*cph+m2*sth*sph+m3*cth
    Omz = 1+n1*sth*cph+n2*sth*sph+n3*cth
    #
    return Omx, Omy, Omz
    
def pathphase (dist, drange, npix, evec, lvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tmesh = np.zeros_like(xmesh)
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
    bevec = evec[3:]
    bevec[2] = bevec[2]*511e3
    bevec = bevec/np.linalg.norm(bevec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        for ycnt in range(ysize):
            opixel = np.array([xmesh[xcnt,ycnt]-x, ymesh[xcnt,ycnt]-y, dist-z])
            opixel = opixel/np.linalg.norm(opixel)
#            print('opixel ')
#            print(np.shape(opixel))
            tmesh[xcnt,ycnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phi = np.arctan2(ymesh - Oy,xmesh - Ox)
    
    #the following only works for electrons with no transverse momentum
    #tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    #phi = np.arctan2(ymesh,xmesh)
    #

#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    #plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #Taylor expansion:
    epsxy = ((xmesh-x)**2 + (ymesh-y)**2)/(dist-z)**2
    if epsxy.all() < 1e-5:
        plmesh = (dist-z)*(epsxy/2)
        print('epsxy: ')
        print(np.min(epsxy))
    else:
        plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    
    incang = np.arccos(np.absolute(np.dot(lvec/np.linalg.norm(lvec),bevec)))

    return tmesh, phi, plmesh, incang
    # returning detector pixel coordinates for trouble shooting/plotting
    #return xmesh, ymesh
    
def pathphase_curved_d (dist, drange, npix, evec, lvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    tmesh = np.zeros_like(xmesh)
    
    #temp angle calculations for z-coordinate of curved detector
    deviation = np.sqrt(xmesh**2 + ymesh**2)
    theta_temp = deviation/dist
    zmesh = dist*np.cos(theta_temp)
    
    
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
    bevec = evec[3:]
    bevec[2] = bevec[2]*511e3
    bevec = bevec/np.linalg.norm(bevec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        for ycnt in range(ysize):
            opixel = np.array([xmesh[xcnt,ycnt]-x, ymesh[xcnt,ycnt]-y, zmesh[xcnt,ycnt]-z])
            opixel = opixel/np.linalg.norm(opixel)
#            print('opixel ')
#            print(np.shape(opixel))
            tmesh[xcnt,ycnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phi = np.arctan2(ymesh - Oy,xmesh - Ox)
    
    #the following only works for electrons with no transverse momentum
    #tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    #phi = np.arctan2(ymesh,xmesh)
    #

#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    #plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #Taylor expansion:
    epsxy = ((xmesh-x)**2 + (ymesh-y)**2)/(dist-z)**2
    if epsxy.all() < 1e-5:
        plmesh = (dist-z)*(epsxy/2)
        print('epsxy: ')
        print(np.min(epsxy))
    else:
        plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    
    incang = np.arccos(np.absolute(np.dot(lvec/np.linalg.norm(lvec),bevec)))

    #return tmesh, phi, plmesh, incang
    # returning detector pixel coordinates for trouble shooting/plotting
    return xmesh, ymesh, zmesh

def pathphase_lin (dist, drange, npix, evec, lvec, polvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    lvec - direction of laser travel
    polvec - normal vector of laser polarization plane
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tmesh = np.zeros_like(xmesh)
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
    bevec = evec[3:]
    bevec[2] = bevec[2]*511e3
    bevec = bevec/np.linalg.norm(bevec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        for ycnt in range(ysize):
            opixel = np.array([xmesh[xcnt,ycnt]-x, ymesh[xcnt,ycnt]-y, dist-z])
            opixel = opixel/np.linalg.norm(opixel)
#            print('opixel ')
#            print(np.shape(opixel))
            tmesh[xcnt,ycnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phi = np.arctan2(ymesh - Oy,xmesh - Ox)
    
    #the following only works for electrons with no transverse momentum
    #tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    #phi = np.arctan2(ymesh,xmesh)
    #

#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    #plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #Taylor expansion:
    epsxy = ((xmesh-x)**2 + (ymesh-y)**2)/(dist-z)**2
    if epsxy.all() < 1e-5:
        plmesh = (dist-z)*(epsxy/2)
        print('epsxy: ')
        print(np.min(epsxy))
    else:
        plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    
#    #angle from laser direction
#    #lvec = np.array([0, 0, 1])
#    incang1 = np.arccos(np.absolute(np.dot(lvec/np.linalg.norm(lvec),bevec)))

#not used:    
    #angle from laser direction, inside plane of polarization
    linplane = np.cross(lvec,polvec)
    incang1 = np.pi/2 - np.arccos(np.absolute(np.dot(linplane/np.linalg.norm(linplane),bevec)))
    
    #angle from laser polarization plane
    #polvec = np.array([1, 0, 0])
    incang2 =np.pi/2 - np.arccos(np.absolute(np.dot(polvec/np.linalg.norm(polvec),bevec)))

    
    return tmesh, phi, plmesh, incang1, incang2
    
    
def pathphase_lin_rotmat (dist, drange, npix, evec, lvec, polvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    lvec - direction of laser travel
    polvec - normal vector of laser polarization plane
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tmesh = np.zeros_like(xmesh)
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
    bevec = evec[3:]
    bevec[2] = bevec[2]*511e3
    bevec = bevec/np.linalg.norm(bevec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        for ycnt in range(ysize):
            opixel = np.array([xmesh[xcnt,ycnt]-x, ymesh[xcnt,ycnt]-y, dist-z])
            opixel = opixel/np.linalg.norm(opixel)
#            print('opixel ')
#            print(np.shape(opixel))
            tmesh[xcnt,ycnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phi = np.arctan2(ymesh - Oy,xmesh - Ox)
    
    #the following only works for electrons with no transverse momentum
    #tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    #phi = np.arctan2(ymesh,xmesh)
    #

#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    #plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #Taylor expansion:
    epsxy = ((xmesh-x)**2 + (ymesh-y)**2)/(dist-z)**2
    if epsxy.all() < 1e-5:
        plmesh = (dist-z)*(epsxy/2)
#        print('epsxy: ')
#        print(np.min(epsxy))
    else:
        plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    ###########################################################################
    ##########   angle calculation via rotation matrix    #####################
    ###########################################################################
    #bevec - beta vector
    
    #zax = np.copy(lvec)
    #zax[2] = -zax[2]
    
    #xax = polvec/np.linalg.norm(polvec)
    
    #temporary: do proper with projections!
    zax = np.array([0.,0.,1.])
    xax = np.array([1.,0.,0.])
    
    bxz = np.copy(bevec)
    bxz[1] = 0
    bxz = bxz/np.linalg.norm(bxz)
#    print('bxz:')
#    print(bxz)
    
    byz = np.copy(bevec)
    byz[0] = 0
    byz = byz/np.linalg.norm(byz)
#    print('byz:')
#    print(byz)

    phi1 = np.arccos(np.abs(np.dot(bxz,zax)))
#    phi1 = phi1[0]
#    print('bxz*zax:')
#    print(np.dot(bxz,zax))
#    print(phi1)
    phi2 = np.arccos(np.abs(np.dot(byz,zax)))
#    phi2 = phi2[0]
#    print('byz*zax:')
#    print(np.dot(byz,zax))
#    print(phi2)
    
    rym = rot_y(phi1)
#    print('roty:')
#    print(rym)
    rxm = rot_x(-phi2)
#    print('rotx:')
#    print(rxm)
    
    rotmat = np.dot(rym,rxm)
#    print('rotmat:')
#    print(rotmat)
        
    xprime = np.dot(rotmat,xax)
#    print('xprime:')
#    print(xprime)
    
#    print('sizes')
#    print(xprime.shape)
#    print(bevec.shape)
#    print(zax.shape)
#    print(xax.shape)
    
#    print(bevec)
#    print(xax)
    
    incang1 = np.arccos(np.abs(np.dot(bevec,zax)))   #theta_I
    incang2 = np.arccos(np.abs(np.dot(xprime,xax)))  #nu
    
#    print('bevec*zax')
#    print(np.dot(bevec,zax))
    
    return tmesh, phi, plmesh, incang1, incang2
    
def pathphase_lin_bvec (dist, drange, npix, evec, lvec, polvec, g0, a0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    lvec - direction of laser travel
    polvec - normal vector of laser polarization plane
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tmesh = np.zeros_like(xmesh)
    #phi = np.zeros_like(xmesh)
    #
    gamma = evec[5]
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
    bevec = evec[3:]
    bevec[2] = bevec[2]*511e3
    nbevec = bevec/np.linalg.norm(bevec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        for ycnt in range(ysize):
            opixel = np.array([xmesh[xcnt,ycnt]-x, ymesh[xcnt,ycnt]-y, dist-z])
            opixel = opixel/np.linalg.norm(opixel)
#            print('opixel ')
#            print(np.shape(opixel))
            tmesh[xcnt,ycnt] = np.arccos(np.absolute(np.dot(opixel,nbevec)))
            #
    #calculating phi
    Ox = x + nbevec[0]*dist
    Oy = y + nbevec[1]*dist
    phi = np.arctan2(ymesh - Oy,xmesh - Ox)
    
    #the following only works for electrons with no transverse momentum
    #tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    #phi = np.arctan2(ymesh,xmesh)
    #

#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    #plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #Taylor expansion:
    epsxy = ((xmesh-x)**2 + (ymesh-y)**2)/(dist-z)**2
    if epsxy.all() < 1e-5:
        plmesh = (dist-z)*(epsxy/2)
#        print('epsxy: ')
#        print(np.min(epsxy))
    else:
        plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    ###########################################################################
    ##########   angle calculation via transversal momentum    #####################
    ###########################################################################
    #bevec - beta vector
    
    h0 = gamma*(1+bevec[2])    
    
    ux0 = gamma*bevec[0]
    uy0 = gamma*bevec[1]
    uzbar = gamma*bevec[2] - a0**2/(4*h0)    
    
    
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uy0**2+uzbar**2))
    nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    
    
    return tmesh, phi, plmesh, thetaI, nu #incang1, incang2


def pathphase_lin_spec (dist, drange, npix, evec, bvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    lvec - direction of laser travel
    polvec - normal vector of laser polarization plane
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(0,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tvec = np.zeros_like(xvec)
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
#    bevec = evec[3:]
#    bevec[2] = bevec[2]*em
#    bevec = bevec/np.linalg.norm(bevec)
    #
    bevec = bvec/np.linalg.norm(bvec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        opixel = np.array([xmesh[xcnt,0]-x, ymesh[xcnt,0]-y, dist-z])
        opixel = opixel/np.linalg.norm(opixel)
        tvec[xcnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phimesh = np.arctan2(ymesh - Oy,xmesh - Ox)
    phi = phimesh[0,npix-1]
        
    return tvec, phi#, plmesh, incang1, incang2

def pathphase_circ_spec (dist, drange, npix, evec, bvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    lvec - direction of laser travel
    polvec - normal vector of laser polarization plane
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    #xvec = np.linspace(0.01,xrange,npix) #lazy workaround to avoid theta = 0
    xvec = np.linspace(0,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tvec = np.zeros_like(xvec)
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
#    bevec = evec[3:]
#    bevec[2] = bevec[2]*em
#    bevec = bevec/np.linalg.norm(bevec)
    #
    bevec = bvec/np.linalg.norm(bvec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        opixel = np.array([xmesh[xcnt,0]-x, ymesh[xcnt,0]-y, dist-z])
        opixel = opixel/np.linalg.norm(opixel)
        tvec[xcnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phimesh = np.arctan2(ymesh - Oy,xmesh - Ox)
    phi = phimesh[0,npix-1]
        
    return tvec, phi#, plmesh, incang1, incang2



def pathphase_pol_grid (theta, phi, npix, evec, lvec, g0):
    '''
    ### UNDER CONSTRUCTION ###
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    x, y, z - position of electron
    xp, yp - transverse momentum of electron
    gamma - lorentz factor of electron
    '''
    #creating detector coordinates:
    xrange = dist*drange/g0
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #zmesh = np.ones_like(xmesh)*dist
    tmesh = np.zeros_like(xmesh)
    #phi = np.zeros_like(xmesh)
    #
    x = evec[0]
    y = evec[1]
    z = evec[2]
    #
    bevec = evec[3:]
    bevec[2] = bevec[2]*511e3
    bevec = bevec/np.linalg.norm(bevec)
#    print('bevec ')
#    print(np.shape(bevec))
    #
    #calculating angles for each pixel
    #including transverse momenta:
    xsize,ysize = np.shape(xmesh)
    for xcnt in range(xsize):
        for ycnt in range(ysize):
            opixel = np.array([xmesh[xcnt,ycnt]-x, ymesh[xcnt,ycnt]-y, dist-z])
            opixel = opixel/np.linalg.norm(opixel)
#            print('opixel ')
#            print(np.shape(opixel))
            tmesh[xcnt,ycnt] = np.arccos(np.absolute(np.dot(opixel,bevec)))
            #
    #calculating phi
    Ox = x + bevec[0]*dist
    Oy = y + bevec[1]*dist
    phi = np.arctan2(ymesh - Oy,xmesh - Ox)
    
    #the following only works for electrons with no transverse momentum
    #tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    #phi = np.arctan2(ymesh,xmesh)
    #

#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    #plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #Taylor expansion:
    epsxy = ((xmesh-x)**2 + (ymesh-y)**2)/(dist-z)**2
    if epsxy.all() < 1e-5:
        plmesh = (dist-z)*(epsxy/2)
        print('epsxy: ')
        print(np.min(epsxy))
    else:
        plmesh = np.sqrt((xmesh-x)**2 + (ymesh-y)**2 + (dist-z)**2)
    
    incang = np.arccos(np.absolute(np.dot(lvec/np.linalg.norm(lvec),bevec)))

    
    return tmesh, phi, plmesh
    

def a0_freq_corr (w0, bet0, a0):
    '''
    correction to base laser frequency to hit resonance with large laser
    intensities
    '''
    wcorr = 4*w0*(1+.5*a0**2)/(1+bet0)**2
    
    return wcorr

def BFF_calc (N0, t_el, warr):
    '''
    N0 - number of electrons
    t_el - time for electron (long. pos./c)
    warr - array with photon frequencies
    '''
#    wmesh, tmesh = np.meshgrid(wvec, tvec)
#    wtmat = np.multiply(wmesh,tmesh)
#das mit meshgrid ist quatsch. es wird teilchenweise bearbeiter -> t ist scalar. w ist schon matrix
    sinsum = 0.
    cossum = 0.
    for cnt in range(N0):
        sinsum += np.sin(warr*t_el[cnt])
        cossum += np.cos(warr*t_el[cnt])
        
    BFF = (1/N0**2)*(sinsum**2 + cossum**2)
    return BFF#np.sqrt(BFF)

    