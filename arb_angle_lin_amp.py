# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 12:25:34 2018

@author: kwl
"""

import numpy as np
from scipy import special as sp

def pathphase (x, y, z, xp, yp):
    '''
    x, y, z - position of electron
    xp, yp - 
    '''

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

#def arb_angle_lin (xvec, yvec, a0, w0, gamma, tinc, nu, N0, nmax, acc):
def arb_angle_lin_amp (evec, tmesh, phi, a0, w0, gamma, bvec, N0, acc, nmax):
#    
#    differential photon distribution
#    linear polarization of incident laser
#    
#    TODO: differentiate between thetaI and nu
#    
#    
#    w - frequency of photons, function parameter
#    thetavec - parametric vector of theta (detector angle)
#    phivec - parametric vector of phi (detector angle)
#    
#    Parameters:
#    a0 - laser parameter
#    N0 - number of laser periods in interaction time
#    T - interaction time
#    w0 - central frequency of incident photons
#    gamma - Lorentz factor of electrons
#    tinc - incidence angle, angle between laser and electron beam
#    theta - polar angle (detector)
#    phi - azimuth angle (detector)
#    nmax - maximum number of harmonics considered
#    ut0 - initial transverse momentum of electron
#    bvec - velocity vector of electron in laser frame, normalized
#    
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
    #rotating beta to interaction angle:
    '''
    rmatz = np.array([[np.cos(nu), -np.sin(nu), 0],
        [np.sin(nu), np.cos(nu), 0],
        [0, 0, 1]])
    #
    rmaty = np.array([[np.cos(tinc), 0, np.sin(tinc)],
        [0, 1, 0],
        [-np.sin(tinc), 0, np.cos(tinc)]])
    #
    #rmat = np.matmul(rmaty, rmatz)
    #rmat = rmaty @ rmatz
    rmat = np.dot(rmaty, rmatz)
    #
    bvec = np.dot(rmat,np.array([[0], [0], [b0]]))
    '''
    bvec = b0*bvec
    print('bvec')
    print(bvec)
    #
    k0 = w0/c
    #k = w/c
    l0 = 2*np.pi*c/w0
    #eta0 = N0*l0/2 #1993
    h0 = gamma*(1+bvec[2])#maybe check again
    r1 = a0/(h0*k0)
    z2 = a0**2/(8*k0*h0**2) #1995
    #z2 = -a0**2/(8*k0*h0**2) #1993
    bzbar = (bvec[2]-(a0**2/(4*gamma*h0)))/(1+(a0**2/(4*gamma*h0)))
    eta0 = (1+bzbar)*N0*l0/2
    #    
    #calculating theta_I:
    ux0 = gamma*bvec[0]
    uy0 = gamma*bvec[1]
    uzbar = gamma*bvec[2] - a0**2/(4*h0)
    #
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uy0**2+uzbar**2))
    nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    #
    print('thetaI:')
    print(thetaI)
    print('nu')
    print(nu)
    #
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    #sti = np.sin(thetaI) #not needed
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    snu = np.sin(nu)
    cnu = np.cos(nu)
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
    ut02 = gamma**2*(bvec[0]**2+bvec[1]**2)
    #
    b1z = (h0**2-(1+.5*a0**2+ut02))/(2*h0**2)
    b1y = bvec[1]/(1+bvec[2])
    b1x = bvec[0]/(1+bvec[2])
    #
    bet1 = np.sqrt(b1x**2+b1y**2+b1z**2)
    #
    Omx, Omy, Omz = co_trafo(thetaI, tmesh, phi, nu)
    #
    g = 1-Omx*b1x-Omy*b1y-Omz*b1z
    #
    Atheta = np.zeros_like(tmesh,dtype = np.complex128)
    Aphi = np.zeros_like(tmesh,dtype = np.complex128)
    #
    #
    for cnt1 in [nmax]:#range(1, nmax+1):
        #resonant w:
        w = cnt1*w0/(1-(1+cth/cti)*b1z)
        #w = w0*2*gamma**2*cnt1*(1+bvec[2])/(1-a0**2/2)#not used here
        k = w/c
        #k = cnt1*k0/(1-b1z-bet1)
        #kbar = k*g-cnt1*k0 + 1e-10
        kbar = 1e-10 #kbar needs to be near zero for resonant w
        #
        #
        #calculate aux functions
        b1 = k*r1*(Omx-Omz*b1x)
        #b1 = cnt1*k0*r1*b1x/bet1 #lin pol extra
        b2 = k*Omz*z2
        #b2 = k*z2*(1+b1z/bet1) #lin pol extra
        #
        #psi0 calculations
        x0 = evec[0]
        y0 = evec[1]
        z0 = evec[2]
        #disabling psi0 for testing purpose
#        print('caution: psi0 disabled!')
#        x0 = 0
#        y0 = 0
#        z0 = 0 
        psi0 = -k*(Omx*x0 + Omy*y0 + Omz*z0)
        #
        B0 = np.zeros_like(tmesh)
        B1 = np.zeros_like(tmesh)
        B2 = np.zeros_like(tmesh)
        for cnt2 in range(-acc, acc+1):
            B0 += 2*sp.jn(cnt2,b2)*sp.jn(cnt1+2*cnt2,b1)
            B1 += (sp.jn(cnt2,b2)*
                (sp.jn(cnt1+2*cnt2+1,b1)+sp.jn(cnt1+2*cnt2-1,b1)))
            B2 += (sp.jn(cnt2,b2)*
                (sp.jn(cnt1+2*cnt2+2,b1)+sp.jn(cnt1+2*cnt2-2,b1)))

        kpref = ec*k/(2*np.pi*np.sqrt(c))
        #integrals as calculated by me:
        #theta coefficients:
        AB0 = (k/k0)*(b1x*(l1*cth*cph+l2*cth*cph-l3*sth)
                + b1y*(m1*cth*cph+m2*cth*sph-m3*sth)
                + b1z*(n1*cth*cph+n2*cth*sph-n3*sth))
        AB1 = k*(r1*(l1*cth*cph+l2*cth*cph-l3*sth)
                -r1*b1x*(n1*cth*cph+n2*cth*sph-n3*sth))
        AB2 = -2*k*(z2*(n1*cth*cph+n2*cth*sph-n3*sth))
        #theta integral:
        Atheta += (np.exp(1j*psi0)*(kpref*(np.sin(kbar*eta0)/kbar)*
                    (AB0*B0 + AB1*B1 + AB2*B2)))
        #
        #phi coefficients:
        BB0 = ((k/k0)*((l2*cph-l1*sph)*b1x + b1y*(m2*cph-m1*sph)
                +b1z*(n2*cph-n1*sph)))
        BB1 = (r1*k*((l2*cph-l1*sph)-b1x*(n2*cph-n1*sph)))
        BB2 = 2*k*z2*(n1*sph-n2*cph)
        #phi integral:
        Aphi += (np.exp(1j*psi0)*(kpref*(np.sin(kbar*eta0)/kbar)*
                    (BB0*B0 + BB1*B1 + BB2*B2)))
        #   
    #
    #tempsum = tempsum1 + tempsum2
    #return (tempsum-np.nanmin(tempsum))/(np.nanmax(tempsum)-np.nanmin(tempsum))
    return Atheta, Aphi, w, psi0
    #    
    #return bvec[0], bvec[1], bvec[2] #output for debugging

def arb_angle_lin_amp2 (evec, bvec, tmesh, phi, a0, w0, gamma, N0, acc, nmax):
#    
#    differential photon distribution
#    linear polarization of incident laser
#    
#    TODO: differentiate between thetaI and nu
#    
#    
#    w - frequency of photons, function parameter
#    thetavec - parametric vector of theta (detector angle)
#    phivec - parametric vector of phi (detector angle)
#    
#    Parameters:
#    a0 - laser parameter
#    N0 - number of laser periods in interaction time
#    T - interaction time
#    w0 - central frequency of incident photons
#    gamma - Lorentz factor of electrons
#    tinc - incidence angle, angle between laser and electron beam
#    theta - polar angle (detector)
#    phi - azimuth angle (detector)
#    nmax - maximum number of harmonics considered
#    ut0 - initial transverse momentum of electron
#    bvec - velocity vector of electron in laser frame, normalized
#    
    #alpha = 0.0072973525664     #fine structure constant
    ec = 1.60217662e-19         #elemental charge
    c = 299792458               #speed of light
    em = 0.51099895e6           #rest mass energy of electreon in eV 
    #
#    print('bvec:')
    #bvec = evec[3:]
    #bvec = np.array(bvec)
#    print(bvec)
    #bvec[2] = bvec[2]*em
#    print(bvec)
    bvec = bvec/c
#    print(bvec)
    bvec = bvec/(np.sqrt((em/c**2)**2 + np.dot(bvec,bvec)/c**2)*c)
    print(bvec)
    #
    k0 = w0/c
    #k = w/c
    l0 = 2*np.pi*c/w0
    #eta0 = N0*l0/2 #1993
    h0 = gamma*(1+bvec[2])#maybe check again
    r1 = a0/(h0*k0)
    z2 = a0**2/(8*k0*h0**2) #1995
    #z2 = -a0**2/(8*k0*h0**2) #1993
    bzbar = (bvec[2]-(a0**2/(4*gamma*h0)))/(1+(a0**2/(4*gamma*h0)))
    eta0 = (1+bzbar)*N0*l0/2
    #    
    #calculating theta_I:
    ux0 = gamma*bvec[0]
    uy0 = gamma*bvec[1]
    uzbar = gamma*bvec[2] - a0**2/(4*h0)
    #
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uy0**2+uzbar**2))
    nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    #
    print('thetaI:')
    print(thetaI)
    print('nu')
    print(nu)
    #
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    #sti = np.sin(thetaI) #not needed
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    snu = np.sin(nu)
    cnu = np.cos(nu)
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
    ut02 = gamma**2*(bvec[0]**2+bvec[1]**2)
#    print('ut02')
#    print(ut02)
    #
    b1z = (h0**2-(1+.5*a0**2+ut02))/(2*h0**2)
    b1y = bvec[1]/(1+bvec[2])
    b1x = bvec[0]/(1+bvec[2])
    #
    bet1 = np.sqrt(b1x**2+b1y**2+b1z**2)
    #
    Omx, Omy, Omz = co_trafo(thetaI, tmesh, phi, nu)
    #
    g = 1-Omx*b1x-Omy*b1y-Omz*b1z
    #
    Atheta = np.zeros_like(tmesh,dtype = np.complex128)
    Aphi = np.zeros_like(tmesh,dtype = np.complex128)
    #
    #
    for cnt1 in [nmax]:#range(1, nmax+1):
        #resonant w:
        #w = cnt1*w0/(1-(1+cth/cti)*b1z)
        #this might cause a problem with theta_I approaching 90 deg
        #
        #this is a simplification that midigates the singularity with the following constraints:
        # gamma0^2 >> a0^2, beta0 ~ 1, theta^2 << 1
        #these should usually be fulfilled
        w = (w0*cnt1*2*(1+bvec[2])*gamma**2)/(1+a0**2/2+(gamma**2)*tmesh**2)
        #
        #w = w0*2*gamma**2*cnt1*(1+bvec[2])/(1-a0**2/2)#not used here
        k = w/c
        #k = cnt1*k0/(1-b1z-bet1)
        #kbar = k*g-cnt1*k0 + 1e-10
        kbar = 1e-10 #kbar needs to be near zero for resonant w
        #
        #
        #calculate aux functions
        b1 = k*r1*(Omx-Omz*b1x)
        #b1 = cnt1*k0*r1*b1x/bet1 #lin pol extra
        b2 = k*Omz*z2
        #b2 = k*z2*(1+b1z/bet1) #lin pol extra
        #
        #psi0 calculations
        x0 = evec[0]
        y0 = evec[1]
        z0 = evec[2]
        #disabling psi0 for testing purpose
#        print('caution: psi0 disabled!')
#        x0 = 0
#        y0 = 0
#        z0 = 0 
        psi0 = -k*(Omx*x0 + Omy*y0 + Omz*z0)
        #
        B0 = np.zeros_like(tmesh)
        B1 = np.zeros_like(tmesh)
        B2 = np.zeros_like(tmesh)
        for cnt2 in range(-acc, acc+1):
            B0 += 2*sp.jn(cnt2,b2)*sp.jn(cnt1+2*cnt2,b1)
            B1 += (sp.jn(cnt2,b2)*
                (sp.jn(cnt1+2*cnt2+1,b1)+sp.jn(cnt1+2*cnt2-1,b1)))
            B2 += (sp.jn(cnt2,b2)*
                (sp.jn(cnt1+2*cnt2+2,b1)+sp.jn(cnt1+2*cnt2-2,b1)))

        #kpref = ec*k/(2*np.pi*np.sqrt(c))
        kpref = ec*k0/(2*np.pi*np.sqrt(c)) #debug test. muss überprüft werden
        #integrals as calculated by me:
        #theta coefficients:
        AB0 = (k/k0)*(b1x*(l1*cth*cph+l2*cth*cph-l3*sth)
                + b1y*(m1*cth*cph+m2*cth*sph-m3*sth)
                + b1z*(n1*cth*cph+n2*cth*sph-n3*sth))
        AB1 = k*(r1*(l1*cth*cph+l2*cth*cph-l3*sth)
                -r1*b1x*(n1*cth*cph+n2*cth*sph-n3*sth))
        AB2 = -2*k*(z2*(n1*cth*cph+n2*cth*sph-n3*sth))
        #theta integral:
        Atheta += (np.exp(1j*psi0)*(kpref*(np.sin(kbar*eta0)/kbar)*
                    (AB0*B0 + AB1*B1 + AB2*B2)))
        #
        #phi coefficients:
        BB0 = ((k/k0)*((l2*cph-l1*sph)*b1x + b1y*(m2*cph-m1*sph)
                +b1z*(n2*cph-n1*sph)))
        BB1 = (r1*k*((l2*cph-l1*sph)-b1x*(n2*cph-n1*sph)))
        BB2 = 2*k*z2*(n1*sph-n2*cph)
        #phi integral:
        Aphi += (np.exp(1j*psi0)*(kpref*(np.sin(kbar*eta0)/kbar)*
                    (BB0*B0 + BB1*B1 + BB2*B2)))
        #   
    #
    #tempsum = tempsum1 + tempsum2
    #return (tempsum-np.nanmin(tempsum))/(np.nanmax(tempsum)-np.nanmin(tempsum))
    return Atheta, Aphi, w, psi0
    
def arb_angle_lin_spec (evec, bvec, wvec, tvec, phi, a0, w0, gamma, N0, acc, nmax):
        #def spec_dist_amp (evec, w, thetavec, phi, a0, w0, gamma, N0, nh, acc):
#    
#    differential photon distribution
#    linear polarization of incident laser
#    
#    
#    
#    w - frequency of photons, function parameter
#    thetavec - parametric vector of theta (detector angle)
#    phivec - parametric vector of phi (detector angle)
#    
#    Parameters:
#    a0 - laser parameter
#    N0 - number of laser periods in interaction time
#    T - interaction time
#    w0 - central frequency of incident photons
#    gamma - Lorentz factor of electrons
#    tinc - incidence angle, angle between laser and electron beam
#    theta - polar angle (detector)
#    phi - azimuth angle (detector)
#    nmax - maximum number of harmonics considered
#    ut0 - initial transverse momentum of electron
#    bvec - velocity vector of electron in laser frame, normalized
#    
    #alpha = 0.0072973525664     #fine structure constant
    ec = 1.60217662e-19         #elemental charge
    c = 299792458               #speed of light
    em = 0.51099895e6           #rest mass energy of electreon in eV 
    #
    #tmesh needs to be created. mesh theta with omega
    w,tmesh = np.meshgrid(wvec,tvec)
    #
    #
    #
    #
    #
#    print('bvec:')
    #bvec = evec[3:]
    bvec = np.array(bvec)
#    print(bvec)
    #bvec[2] = bvec[2]*em
#    print(bvec)
    bvec = bvec/c
#    print(bvec)
    bvec = bvec/(np.sqrt((em/c**2)**2 + np.dot(bvec,bvec)/c**2)*c)
#    print(bvec)
    #
    k = w/c
    k0 = w0/c
    #k = w/c
    l0 = 2*np.pi*c/w0
    #eta0 = N0*l0/2 #1993
    h0 = gamma*(1+bvec[2])#maybe check again
    r1 = a0/(h0*k0)
    z2 = a0**2/(8*k0*h0**2) #1995
    #z2 = -a0**2/(8*k0*h0**2) #1993
    bzbar = (bvec[2]-(a0**2/(4*gamma*h0)))/(1+(a0**2/(4*gamma*h0)))
    eta0 = (1+bzbar)*N0*l0/2
    #eta0 = N0*l0 #1993
    #    
    #calculating theta_I:
    ux0 = gamma*bvec[0]
    uy0 = gamma*bvec[1]
    uzbar = gamma*bvec[2] - a0**2/(4*h0)
    #
    thetaI = np.arccos(uzbar/np.sqrt(ux0**2+uy0**2+uzbar**2))
    nu = np.arccos(uzbar/np.sqrt(ux0**2+uzbar**2))
    #
    print('thetaI:')
    print(thetaI)
    print('nu')
    print(nu)
    #
    cti = np.cos(thetaI)
    cth = np.cos(tmesh)
    cph = np.cos(phi)
    #sti = np.sin(thetaI) #not needed
    sth = np.sin(tmesh)
    sph = np.sin(phi)
    snu = np.sin(nu)
    cnu = np.cos(nu)
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
    ut02 = gamma**2*(bvec[0]**2+bvec[1]**2)
#    print('ut02')
#    print(ut02)
    #
    b1z = (h0**2-(1+.5*a0**2+ut02))/(2*h0**2)
    b1y = bvec[1]/(1+bvec[2])
    b1x = bvec[0]/(1+bvec[2])
    #
    #
    bet1 = np.sqrt(b1x**2+b1y**2+b1z**2)
    #
    Omx, Omy, Omz = co_trafo(thetaI, tmesh, phi, nu)
    #
    g = 1-Omx*b1x-Omy*b1y-Omz*b1z
    #
    Atheta = np.zeros_like(tmesh,dtype = np.complex128)
    Aphi = np.zeros_like(tmesh,dtype = np.complex128)
    #
    #
    for cnt1 in range(1, nmax+1):#[nmax]:
        #resonant w:
        #w = cnt1*w0/(1-(1+cth/cti)*b1z)
        #w = w0*2*gamma**2*cnt1*(1+bvec[2])/(1-a0**2/2)#not used here
        #k = cnt1*k0/(1-b1z-bet1)
        kbar = k*g-cnt1*k0
        #kbar = 1e-10 #kbar needs to be near zero for resonant w
        #
        #
        #calculate aux functions
        b1 = k*r1*(Omx-Omz*b1x)
        #b1 = cnt1*k0*r1*b1x/bet1 #lin pol extra
        b2 = k*Omz*z2
        #b2 = k*z2*(1+b1z/bet1) #lin pol extra
        #
        #psi0 calculations
        x0 = evec[0]
        y0 = evec[1]
        z0 = evec[2]
        #disabling psi0 for testing purpose
#        print('caution: psi0 disabled!')
#        x0 = 0
#        y0 = 0
#        z0 = 0 
        psi0 = -k*(Omx*x0 + Omy*y0 + Omz*z0)
        #
        B0 = np.zeros_like(tmesh)
        B1 = np.zeros_like(tmesh)
        B2 = np.zeros_like(tmesh)
        for cnt2 in range(-acc, acc+1):
            B0 += 2*sp.jn(cnt2,b2)*sp.jn(cnt1+2*cnt2,b1)
            B1 += (sp.jn(cnt2,b2)*
                (sp.jn(cnt1+2*cnt2+1,b1)+sp.jn(cnt1+2*cnt2-1,b1)))
            B2 += (sp.jn(cnt2,b2)*
                (sp.jn(cnt1+2*cnt2+2,b1)+sp.jn(cnt1+2*cnt2-2,b1)))

        #kpref = ec*k/(2*np.pi*np.sqrt(c))
        kpref = ec*k0/(2*np.pi*np.sqrt(c))
        #integrals as calculated by me:
        #theta coefficients:
        
        print('b1x:')
        print(b1x)
        print('b1y:')
        print(b1y)
        print('b1z:')
        print(b1z)
        AB0 = (k/k0)*(b1x*(l1*cth*cph+l2*cth*cph-l3*sth)
                + b1y*(m1*cth*cph+m2*cth*sph-m3*sth)
                + b1z*(n1*cth*cph+n2*cth*sph-n3*sth))
        AB1 = k*(r1*(l1*cth*cph+l2*cth*cph-l3*sth)
                -r1*b1x*(n1*cth*cph+n2*cth*sph-n3*sth))
        AB2 = -2*k*(z2*(n1*cth*cph+n2*cth*sph-n3*sth))
        #
#        AB0 = 0
#        print('theta B0 coeff:')
#        print(AB0)
        #theta integral:
        Atheta += (np.exp(1j*psi0)*(kpref*(np.sin(kbar*eta0)/kbar)*
                    (AB0*B0 + AB1*B1 + AB2*B2)))
        #
        #phi coefficients:
        BB0 = ((k/k0)*((l2*cph-l1*sph)*b1x + b1y*(m2*cph-m1*sph)
                +b1z*(n2*cph-n1*sph)))
        BB1 = (r1*k*((l2*cph-l1*sph)-b1x*(n2*cph-n1*sph)))
        BB2 = 2*k*z2*(n1*sph-n2*cph)
        #
#        print('phi B0 coeff:')
#        print(BB0)
        #phi integral:
        Aphi += (np.exp(1j*psi0)*(kpref*(np.sin(kbar*eta0)/kbar)*
                    (BB0*B0 + BB1*B1 + BB2*B2)))
        #   
    #
    #tempsum = tempsum1 + tempsum2
    #return (tempsum-np.nanmin(tempsum))/(np.nanmax(tempsum)-np.nanmin(tempsum))
    return Atheta, Aphi#, w, psi0