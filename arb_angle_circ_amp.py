# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 11:25:42 2018

@author: kwl


eliminating infinite sums with Graf's theorem
"""

import numpy as np
from scipy import special as sp

def pathphase (dist, drange, npix, x, y, z, xp, yp, gamma, g0):
    '''
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
    xrange = g0*drange/dist
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #
    #calculating angles for each pixel
    tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    phi = np.arctan2(ymesh,xmesh)
    #
#    print('xmax:')
#    print(np.max(xmesh/g0))
    #calculating distance to each pixel
    plmesh = np.sqrt((xmesh/g0-x)**2 + (ymesh/g0-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    #taylor series: (temp fix for only z)
#    xyterms = (xmesh/g0+x)**2 + (ymesh/g0+y)**2
#    plmesh = dist*z/np.sqrt(dist**2 + xyterms)
#    print(dist**2 + np.max(xyterms))
#    print(dist**2 + np.min(xyterms))
#    print(np.max(plmesh))
#    print(np.min(plmesh))
    return tmesh, phi, plmesh
    
def pathphase2 (dist, drange, npix, earray, Ne):
    '''
    this one handles the whole list of electrons
    reference particle must be first in the list    
    
    dist -  distance from nominal particle interaction to center of detector
            (must remain constant for every electron in bunch)
    drange -    size of the detector: square detector
                range of angle gamma*x/z and gamma*y/z: -drange to +drange
    npix -  number of pixels in the detector:
            square detector: total number of pixels: npix*npix
    earray - list of electron coordinates and momenta
    Ne - number of considered electrons
    '''
    g0 = earray[0, 5]
    #creating detector coordinates:
    xrange = g0*drange/dist
    xvec = np.linspace(-xrange,xrange,npix)
    #
    xmesh,ymesh = np.meshgrid(xvec,xvec)
    #
    #calculating angles for each pixel
    tmesh = np.sqrt((dist*(xmesh/g0-x))**2+(dist*(ymesh/g0-y))**2)/(dist+z)
    phi = np.arctan2(ymesh,xmesh)
    #
    #calculating distance to each pixel
    plmesh = np.sqrt((xmesh/g0-x)**2 + (ymesh/g0-y)**2 + (dist-z)**2)
    #this might need to be Taylored
    return tmesh, phi, plmesh
    

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

#def arb_angle_circ_amp (x0, y0, z0, tmesh, phi, a0, w0, gamma, incangle, N0, nmax):#, acc):
def arb_angle_circ_amp (evec, tmesh, phi, a0, w0, gamma, incangle, N0, nmax):#, acc):
    '''
    differential photon distribution
    circular polarization of incident laser
    
    evec - 6-vector: coordinate and momenta of electron
    
    w - frequency of photons, function parameter
    thetavec - parametric vector of theta (detector angle)
    phivec - parametric vector of phi (detector angle)
    
    Parameters:
    a0 - laser parameter
    N0 - number of laser periods in interaction time
    T - interaction time
    w0 - central frequency of incident photons
    gamma - Lorentz factor of electrons
    thetaI - incidence angle, angle between laser and electron beam
    theta - polar angle (detector)
    phi - azimuth angle (detector)
    nmax - maximum number of harmonics considered
    ut0 - initial transverse momentum of electron
    bvec - velocity vector of electron in laser frame
    '''
    #alpha = 0.0072973525664     #fine structure constant
    ec = 1.60217662e-19         #elemental charge
    c = 299792458               #speed of light
    #
#    xmesh,ymesh = np.meshgrid(xvec,yvec)
#    #
#    tmesh = np.sqrt((xmesh/gamma)**2+(ymesh/gamma)**2)
#    phi = np.arctan2(ymesh,xmesh)
    #
    b0 = np.sqrt(1-1/gamma**2)
    #
    '''
    TODO: consider transverse momentum in incangle
    '''
    #rotating beta to interaction angle:
    bvec = b0*np.array([[np.sin(incangle)], [0], [np.cos(incangle)]])
    #
    k0 = w0/c
    #k = w/c
    l0 = 2*np.pi*c/w0
    #eta0 = N0*l0
    h0 = gamma*(1+bvec[2])#maybe check again
    r1 = a0/(h0*k0)
    #
    bzbar = (bvec[2]-(a0**2/(4*gamma*h0)))/(1+(a0**2/(4*gamma*h0)))
    eta0 = (1+bzbar)*N0*l0/2
    #    
    #calculating theta_I:
    #for circular polarization:
    #beta_y = 0
    #therefore thetaI = nu
    #
    ux0 = gamma*bvec[0]
    #uy0 = gamma*bvec[1]
    uzbar = gamma*bvec[2] - a0**2/(4*h0)
    #
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    #nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    nu = thetaI    
    #
    #'''
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    sti = np.sin(thetaI)
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    #for circ. pol. thetaI = nu
    cnu = cti
    snu = sti
    #'''
    #
    l1 = cnu
    l2 = -snu*np.sqrt(1-cti**2/cnu**2)
    l3 = np.tan(nu)*cti
    #m not used
    m1 = 0
    m2 = cti/cnu
    m3 = np.sqrt(1-cti**2/cnu**2)
    n1 = -snu
    n2 = -cnu*np.sqrt(1-cti**2/cnu**2)
    n3 = cti
    #
    Omx, Omy, Omz = co_trafo(thetaI, tmesh, phi, thetaI)
    #
    ut02 = gamma**2*(bvec[0]**2+bvec[1]**2)
    #
    b1x = bvec[0]/(1+bvec[2])
    b1y = bvec[1]/(1+bvec[2])
    b1z = (h0**2-(1+.5*a0**2+ut02))/(2*h0**2)
    #g = 1-(1+cth/cti)*b1z #maybe wrong
    #g = 1-b1x*(cti*sth*cph+sti*cth)-b1z*(1-sti*sth*cph+cti*cth)
    g = 1-Omx*b1x-Omy*b1y-Omz*b1z
    #
    #
    #
#    tempsum1 = np.zeros_like(tmesh)
#    tempsum2 = np.zeros_like(tmesh)
    Atheta = np.zeros_like(tmesh,dtype = np.complex128)
    Aphi = np.zeros_like(tmesh,dtype = np.complex128)
    I0 = np.zeros_like(tmesh,dtype = np.complex128)
    Is = np.zeros_like(tmesh,dtype = np.complex128)
    Ic = np.zeros_like(tmesh,dtype = np.complex128)
    for cnt1 in [nmax]:#range(1, nmax+1):
        #w = w0*2*gamma**2*cnt1*(1+bvec[2])/(1-a0**2/2)
        w = cnt1*w0/(1-(1+cth/cti)*b1z) #omega as resonant function
        k = w/c
        #kbar = k*g-cnt1*k0
        kbar = 1e-5 #kbar needs to be near zero for resonant w
        #
        #bs = ((k*r1/np.sqrt(2))*(cth*sti*(1-b1z)-b1x)
        #    +sth*cph*(cti+b1x*sti)) #maybe wrong
        #bs = ((k*r1/np.sqrt(2))*(cth*sti+sth*cph*(cti+b1x*sti)-b1x
        #    -b1x*cti*cth))
        bs = k*r1*(Omx-Omz*b1x)/np.sqrt(2)
        #bc = (k*r1/np.sqrt(2))*sth*cph
        bc = k*r1*(Omy-Omz*b1y)/np.sqrt(2)
        bt2 = bs**2+bc**2
        bt = np.sqrt(bt2)
        #
        x0 = evec[0]
        y0 = evec[1]
        z0 = evec[2]
        #disabling psi0 for testing purpose
#        print('caution: psi0 disabled!')
#        x0 = 0
#        y0 = 0
#        z0 = 0 
        psi0 = -k*(Omx*x0 + Omy*y0 + Omz*z0)
        #xtheta = cnu*cth*cph - snu*sth
        xtheta = l1*cth*cph +l2*cth*sph -l3*sth
        #ytheta = cth*sph
        ytheta = m1*cth*cph + m2*cth*sph -m3*sth
        #ztheta = cnu*cth*cph - snu*sth
        ztheta = n1*cth*cph +n2*cth*sph - n3*sth
        I0t = b1x*xtheta + b1y*ytheta + b1z*ztheta
        Ict = (k0*r1/np.sqrt(2))*(xtheta - b1x*ztheta)
        Ist = (k0*r1/np.sqrt(2))*(ytheta - b1y*ztheta)
        #
        xphi = -cnu*sph
        yphi = cph
        zphi = snu*sph
        I0p = b1x*xphi + b1y*yphi + b1z*zphi
        Icp = (k0*r1/np.sqrt(2))*(xphi - b1x*zphi)
        Isp = (k0*r1/np.sqrt(2))*(yphi - b1y*zphi)
        #
        I0 += np.exp(1j*psi0)*2*sp.jn(cnt1,bt)*((bs + 1j*bc)/bt)**cnt1

        Ic += np.exp(1j*psi0)*(sp.jn(cnt1+1,bt)*((bs + 1j*bc)/bt)**(cnt1+1) +
                sp.jn(cnt1-1,bt)*((bs + 1j*bc)/bt)**(cnt1-1))

        Is += -1j*np.exp(1j*psi0)*(sp.jn(cnt1+1,bt)*((bs + 1j*bc)/bt)**(cnt1+1) -
                sp.jn(cnt1-1,bt)*((bs + 1j*bc)/bt)**(cnt1-1))
                   
#        tempsum1 += ((ec**2*k**2/(4*np.pi**2*c))*(np.sin(kbar*eta0)/kbar)**2*
#                    np.absolute((I0t*I0 + Ict*Ic + Ist*Is))**2)
#        tempsum2 += ((ec**2*k**2/(4*np.pi**2*c))*(np.sin(kbar*eta0)/kbar)**2*
#                    np.absolute((I0p*I0 +Icp*Ic + Isp*Is))**2)
        kpref = ec*k/(2*np.pi*np.sqrt(c))
        Atheta += kpref*(I0t*I0 + Ict*Ic + Ist*Is)*(np.sin(kbar*eta0)/kbar)
        Aphi += kpref*(I0p*I0 +Icp*Ic + Isp*Is)*(np.sin(kbar*eta0)/kbar)

        
    #
    #    tempsum =  tempsum1 + tempsum2
    #return (tempsum-np.nanmin(tempsum))/(np.nanmax(tempsum)-np.nanmin(tempsum))
    return Atheta, Aphi, w, psi0
    #return bvec #output for debugging
    
def arb_angle_circ_amp_nophase (evec, tmesh, phi, a0, w0, gamma, incangle, N0, nmax):#, acc):
    '''
    differential photon distribution
    circular polarization of incident laser
    
    evec - 6-vector: coordinate and momenta of electron
    
    w - frequency of photons, function parameter
    thetavec - parametric vector of theta (detector angle)
    phivec - parametric vector of phi (detector angle)
    
    Parameters:
    a0 - laser parameter
    N0 - number of laser periods in interaction time
    T - interaction time
    w0 - central frequency of incident photons
    gamma - Lorentz factor of electrons
    thetaI - incidence angle, angle between laser and electron beam
    theta - polar angle (detector)
    phi - azimuth angle (detector)
    nmax - maximum number of harmonics considered
    ut0 - initial transverse momentum of electron
    bvec - velocity vector of electron in laser frame
    '''
    #alpha = 0.0072973525664     #fine structure constant
    ec = 1.60217662e-19         #elemental charge
    c = 299792458               #speed of light
    #
#    xmesh,ymesh = np.meshgrid(xvec,yvec)
#    #
#    tmesh = np.sqrt((xmesh/gamma)**2+(ymesh/gamma)**2)
#    phi = np.arctan2(ymesh,xmesh)
    #
    b0 = np.sqrt(1-1/gamma**2)
    #
    '''
    TODO: consider transverse momentum in incangle
    '''
    #rotating beta to interaction angle:
    bvec = b0*np.array([[np.sin(incangle)], [0], [np.cos(incangle)]])
    #
    k0 = w0/c
    #k = w/c
    l0 = 2*np.pi*c/w0
    #eta0 = N0*l0
    h0 = gamma*(1+bvec[2])#maybe check again
    r1 = a0/(h0*k0)
    #
    bzbar = (bvec[2]-(a0**2/(4*gamma*h0)))/(1+(a0**2/(4*gamma*h0)))
    eta0 = (1+bzbar)*N0*l0/2
    #    
    #calculating theta_I:
    #for circular polarization:
    #beta_y = 0
    #therefore thetaI = nu
    #
    ux0 = gamma*bvec[0]
    #uy0 = gamma*bvec[1]
    uzbar = gamma*bvec[2] - a0**2/(4*h0)
    #
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    #nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    nu = thetaI    
    #
    #'''
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    sti = np.sin(thetaI)
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    #for circ. pol. thetaI = nu
    cnu = cti
    snu = sti
    #'''
    #
    l1 = cnu
    l2 = -snu*np.sqrt(1-cti**2/cnu**2)
    l3 = np.tan(nu)*cti
    #m not used
    m1 = 0
    m2 = cti/cnu
    m3 = np.sqrt(1-cti**2/cnu**2)
    n1 = -snu
    n2 = -cnu*np.sqrt(1-cti**2/cnu**2)
    n3 = cti
    #
    Omx, Omy, Omz = co_trafo(thetaI, tmesh, phi, thetaI)
    #
    ut02 = gamma**2*(bvec[0]**2+bvec[1]**2)
    #
    b1x = bvec[0]/(1+bvec[2])
    b1y = bvec[1]/(1+bvec[2])
    b1z = (h0**2-(1+.5*a0**2+ut02))/(2*h0**2)
    #g = 1-(1+cth/cti)*b1z #maybe wrong
    #g = 1-b1x*(cti*sth*cph+sti*cth)-b1z*(1-sti*sth*cph+cti*cth)
    g = 1-Omx*b1x-Omy*b1y-Omz*b1z
    #
    #
    #
#    tempsum1 = np.zeros_like(tmesh)
#    tempsum2 = np.zeros_like(tmesh)
    Atheta = np.zeros_like(tmesh,dtype = np.complex128)
    Aphi = np.zeros_like(tmesh,dtype = np.complex128)
    I0 = np.zeros_like(tmesh,dtype = np.complex128)
    Is = np.zeros_like(tmesh,dtype = np.complex128)
    Ic = np.zeros_like(tmesh,dtype = np.complex128)
    for cnt1 in [nmax]:#range(1, nmax+1):
        #w = w0*2*gamma**2*cnt1*(1+bvec[2])/(1-a0**2/2)
        w = cnt1*w0/(1-(1+cth/cti)*b1z) #omega as resonant function
        k = w/c
        #kbar = k*g-cnt1*k0
        kbar = 1e-5 #kbar needs to be near zero for resonant w
        #
        #bs = ((k*r1/np.sqrt(2))*(cth*sti*(1-b1z)-b1x)
        #    +sth*cph*(cti+b1x*sti)) #maybe wrong
        #bs = ((k*r1/np.sqrt(2))*(cth*sti+sth*cph*(cti+b1x*sti)-b1x
        #    -b1x*cti*cth))
        bs = k*r1*(Omx-Omz*b1x)/np.sqrt(2)
        #bc = (k*r1/np.sqrt(2))*sth*cph
        bc = k*r1*(Omy-Omz*b1y)/np.sqrt(2)
        bt2 = bs**2+bc**2
        bt = np.sqrt(bt2)
        #
#        x0 = evec[0]
#        y0 = evec[1]
#        z0 = evec[2]
        #disabling psi0 for testing purpose
        print('caution: psi0 disabled!')
        x0 = 0
        y0 = 0
        z0 = 0 
        psi0 = -k*(Omx*x0 + Omy*y0 + Omz*z0)
        #xtheta = cnu*cth*cph - snu*sth
        xtheta = l1*cth*cph +l2*cth*sph -l3*sth
        #ytheta = cth*sph
        ytheta = m1*cth*cph + m2*cth*sph -m3*sth
        #ztheta = cnu*cth*cph - snu*sth
        ztheta = n1*cth*cph +n2*cth*sph - n3*sth
        I0t = b1x*xtheta + b1y*ytheta + b1z*ztheta
        Ict = (k0*r1/np.sqrt(2))*(xtheta - b1x*ztheta)
        Ist = (k0*r1/np.sqrt(2))*(ytheta - b1y*ztheta)
        #
        xphi = -cnu*sph
        yphi = cph
        zphi = snu*sph
        I0p = b1x*xphi + b1y*yphi + b1z*zphi
        Icp = (k0*r1/np.sqrt(2))*(xphi - b1x*zphi)
        Isp = (k0*r1/np.sqrt(2))*(yphi - b1y*zphi)
        #
        I0 += np.exp(1j*psi0)*2*sp.jn(cnt1,bt)*((bs + 1j*bc)/bt)**cnt1

        Ic += np.exp(1j*psi0)*(sp.jn(cnt1+1,bt)*((bs + 1j*bc)/bt)**(cnt1+1) +
                sp.jn(cnt1-1,bt)*((bs + 1j*bc)/bt)**(cnt1-1))

        Is += -1j*np.exp(1j*psi0)*(sp.jn(cnt1+1,bt)*((bs + 1j*bc)/bt)**(cnt1+1) -
                sp.jn(cnt1-1,bt)*((bs + 1j*bc)/bt)**(cnt1-1))
                   
#        tempsum1 += ((ec**2*k**2/(4*np.pi**2*c))*(np.sin(kbar*eta0)/kbar)**2*
#                    np.absolute((I0t*I0 + Ict*Ic + Ist*Is))**2)
#        tempsum2 += ((ec**2*k**2/(4*np.pi**2*c))*(np.sin(kbar*eta0)/kbar)**2*
#                    np.absolute((I0p*I0 +Icp*Ic + Isp*Is))**2)
        kpref = ec*k/(2*np.pi*np.sqrt(c))
        Atheta += kpref*(I0t*I0 + Ict*Ic + Ist*Is)*(np.sin(kbar*eta0)/kbar)
        Aphi += kpref*(I0p*I0 +Icp*Ic + Isp*Is)*(np.sin(kbar*eta0)/kbar)

        
    #
    #    tempsum =  tempsum1 + tempsum2
    #return (tempsum-np.nanmin(tempsum))/(np.nanmax(tempsum)-np.nanmin(tempsum))
    return Atheta, Aphi, w, psi0
    #return bvec #output for debugging
    
    
def arb_angle_circ_spec (evec, bvec, wvec, tvec, phi, a0, w0, gamma, N0, nmax):
    '''
    differential photon distribution
    circular polarization of incident laser
    
    evec - 6-vector: coordinate and momenta of electron
    
    w - frequency of photons, function parameter
    thetavec - parametric vector of theta (detector angle)
    phivec - parametric vector of phi (detector angle)
    
    Parameters:
    a0 - laser parameter
    N0 - number of laser periods in interaction time
    T - interaction time
    w0 - central frequency of incident photons
    gamma - Lorentz factor of electrons
    thetaI - incidence angle, angle between laser and electron beam
    theta - polar angle (detector)
    phi - azimuth angle (detector)
    nmax - maximum number of harmonics considered
    ut0 - initial transverse momentum of electron
    bvec - velocity vector of electron in laser frame
    '''
    #alpha = 0.0072973525664     #fine structure constant
    ec = 1.60217662e-19         #elemental charge
    c = 299792458               #speed of light
    em = 0.51099895e6           #rest mass energy of electreon in eV 
    #
#    xmesh,ymesh = np.meshgrid(xvec,yvec)
#    #
#    tmesh = np.sqrt((xmesh/gamma)**2+(ymesh/gamma)**2)
#    phi = np.arctan2(ymesh,xmesh)
    #
    w,tmesh = np.meshgrid(wvec,tvec)
    #
    #
    #
    bvec = bvec/c
    bvec = bvec/(np.sqrt((em/c**2)**2 + np.dot(bvec,bvec)/c**2)*c)
    #
    '''
    TODO: consider transverse momentum in incangle
    '''
    #rotating beta to interaction angle:
    #
    k0 = w0/c
    #k = w/c
    l0 = 2*np.pi*c/w0
    #eta0 = N0*l0
    h0 = gamma*(1+bvec[2])#maybe check again
    r1 = a0/(h0*k0)
    #
    bzbar = (bvec[2]-(a0**2/(4*gamma*h0)))/(1+(a0**2/(4*gamma*h0)))
    eta0 = (1+bzbar)*N0*l0/2
    #    
    #calculating theta_I:
    #for circular polarization:
    #beta_y = 0
    #therefore thetaI = nu
    #
    ux0 = gamma*bvec[0]
    #uy0 = gamma*bvec[1]
    uzbar = gamma*bvec[2] - a0**2/(4*h0)
    #
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    #nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    nu = thetaI    
    #
    #'''
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    sti = np.sin(thetaI)
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    #for circ. pol. thetaI = nu
    cnu = cti
    snu = sti
    #'''
    #
    l1 = cnu
    l2 = -snu*np.sqrt(1-cti**2/cnu**2)
    l3 = np.tan(nu)*cti
    #m not used
    m1 = 0
    m2 = cti/cnu
    m3 = np.sqrt(1-cti**2/cnu**2)
    n1 = -snu
    n2 = -cnu*np.sqrt(1-cti**2/cnu**2)
    n3 = cti
    #
    Omx, Omy, Omz = co_trafo(thetaI, tmesh, phi, thetaI)
    #
    ut02 = gamma**2*(bvec[0]**2+bvec[1]**2)
    #
    b1x = bvec[0]/(1+bvec[2])
    b1y = bvec[1]/(1+bvec[2])
    b1z = (h0**2-(1+.5*a0**2+ut02))/(2*h0**2)
    #g = 1-(1+cth/cti)*b1z #maybe wrong
    #g = 1-b1x*(cti*sth*cph+sti*cth)-b1z*(1-sti*sth*cph+cti*cth)
    g = 1-Omx*b1x-Omy*b1y-Omz*b1z
    #
    #
    #
#    tempsum1 = np.zeros_like(tmesh)
#    tempsum2 = np.zeros_like(tmesh)
    Atheta = np.zeros_like(tmesh,dtype = np.complex128)
    Aphi = np.zeros_like(tmesh,dtype = np.complex128)
    I0 = np.zeros_like(tmesh,dtype = np.complex128)
    Is = np.zeros_like(tmesh,dtype = np.complex128)
    Ic = np.zeros_like(tmesh,dtype = np.complex128)
    intpart = np.zeros_like(tmesh)
    for cnt1 in [nmax]:#range(1, nmax+1):#
        #w = w0*2*gamma**2*cnt1*(1+bvec[2])/(1-a0**2/2)
        #w = cnt1*w0/(1-(1+cth/cti)*b1z) #omega as resonant function
        k = w/c
        kbar = k*g-cnt1*k0
        #kbar = 1e-5 #kbar needs to be near zero for resonant w
        #
        #bs = ((k*r1/np.sqrt(2))*(cth*sti*(1-b1z)-b1x)
        #    +sth*cph*(cti+b1x*sti)) #maybe wrong
        #bs = ((k*r1/np.sqrt(2))*(cth*sti+sth*cph*(cti+b1x*sti)-b1x
        #    -b1x*cti*cth))
        bs = k*r1*(Omx-Omz*b1x)/np.sqrt(2)
        #bc = (k*r1/np.sqrt(2))*sth*cph
        bc = k*r1*(Omy-Omz*b1y)/np.sqrt(2)
        bt2 = bs**2+bc**2
        bt = np.sqrt(bt2)
        #
        x0 = evec[0]
        y0 = evec[1]
        z0 = evec[2]
        #disabling psi0 for testing purpose
#        print('caution: psi0 disabled!')
#        x0 = 0
#        y0 = 0
#        z0 = 0 
        psi0 = -k*(Omx*x0 + Omy*y0 + Omz*z0)
        #xtheta = cnu*cth*cph - snu*sth
        xtheta = l1*cth*cph +l2*cth*sph -l3*sth
        #ytheta = cth*sph
        ytheta = m1*cth*cph + m2*cth*sph -m3*sth
        #ztheta = cnu*cth*cph - snu*sth
        ztheta = n1*cth*cph +n2*cth*sph - n3*sth
        I0t = b1x*xtheta + b1y*ytheta + b1z*ztheta
        Ict = (k0*r1/np.sqrt(2))*(xtheta - b1x*ztheta)
        Ist = (k0*r1/np.sqrt(2))*(ytheta - b1y*ztheta)
        #
        xphi = -cnu*sph
        yphi = cph
        zphi = snu*sph
        I0p = b1x*xphi + b1y*yphi + b1z*zphi
        Icp = (k0*r1/np.sqrt(2))*(xphi - b1x*zphi)
        Isp = (k0*r1/np.sqrt(2))*(yphi - b1y*zphi)
        #
        #print(bt)
        I0 += np.exp(1j*psi0)*2*sp.jn(cnt1,bt)*((bs + 1j*bc)/bt)**cnt1

        Ic += np.exp(1j*psi0)*(sp.jn(cnt1+1,bt)*((bs + 1j*bc)/bt)**(cnt1+1) +
                sp.jn(cnt1-1,bt)*((bs + 1j*bc)/bt)**(cnt1-1))

        Is += -1j*np.exp(1j*psi0)*(sp.jn(cnt1+1,bt)*((bs + 1j*bc)/bt)**(cnt1+1) -
                sp.jn(cnt1-1,bt)*((bs + 1j*bc)/bt)**(cnt1-1))
                   
        kpref = ec*k/(2*np.pi*np.sqrt(c))
        Atheta += kpref*(I0t*I0 + Ict*Ic + Ist*Is)*(np.sin(kbar*eta0)/kbar)
        Aphi += kpref*(I0p*I0 +Icp*Ic + Isp*Is)*(np.sin(kbar*eta0)/kbar)
        
        intpart += np.absolute(Atheta)**2 + np.absolute(Aphi)**2
        #######################################################################
        ##### theta = 0 case       ############################################
        #######################################################################
        t0index = tmesh == 0
        if np.any(t0index):
            if thetaI == 0:
                intp = ((ec*k0)**2/(np.pi**2*c)*(np.sin(kbar*eta0)/kbar)**2)/2
                intpart[t0index] = intp[t0index]
            else:
                bs = k*r1*(sti-(1+cti)*b1x)/np.sqrt(2)
                Jprime = (sp.jn(cnt1-1,bs)-sp.jn(cnt1+1,bs))/2
                intp = (kpref**2*(np.sin(kbar*eta0)/kbar)**2
                    *(a0**2/(2*h0**2))*(Jprime**2 + (cti+b1x*sti)**2
                    *(cnt1*sp.jn(cnt1,bs)/bs)**2))
                intpart[t0index] = intp[t0index]
    #
    #    tempsum =  tempsum1 + tempsum2
    #return (tempsum-np.nanmin(tempsum))/(np.nanmax(tempsum)-np.nanmin(tempsum))
    #return Atheta, Aphi#, w, psi0
    return intpart
    
    