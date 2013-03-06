#!/usr/bin/env python
# -*- coding: ascii -*- 
""" AstrePhase - Module d'affichage de la phase d'un astre vue depuis la Terre a une date donnee
Modules utilises : python (2.7.3), pyephem (3.7.5.1), matplotlib (1.1.1), numpy (1.6.2)
"""
__author__ = "J-P Fontaine & J-B Lekien"
__copyright__ = "Copyright (c) 2012, The Practice Astronom Guide Book Project"
__date__ = "2012/10/20"
__version__ = "1.0.0"

import ephem
import pylab
import math
import numpy
import matplotlib.patches as patches

def compute_khi(obs, ra_dec, astre) :
    """ Calcul du khi representant l'angle d'inclinaison du terminateur
    ra_dec = True  : Astrometric geocentric position represented by the properties : a_ra and a_dec.
    ra_dec = ...   : Apparent geocentric position represented by the properties    : g_ra and g_dec
    ra_dec = False : Apparent topocentric position represented by the properties   : ra and dec, and alt and az.
    """
    soleil = ephem.Sun()
    if soleil.name == astre.name :
        return 0.
    astre.compute(obs)
    soleil.compute(obs)
    if ra_dec :
        numerateur = math.cos(soleil.a_dec) * math.sin(soleil.a_ra - astre.a_ra)
        denominateur = math.cos(astre.a_dec) * math.sin(soleil.a_dec) - math.sin(astre.a_dec) * math.cos(soleil.a_dec) * math.cos(soleil.a_ra - astre.a_ra)
    else :
        numerateur = - math.cos(soleil.alt) * math.sin(soleil.az - astre.az)
        denominateur = math.cos(astre.alt) * math.sin(soleil.alt) - math.sin(astre.alt) * math.cos(soleil.alt) * math.cos(soleil.az - astre.az)
    khi = numpy.arctan( numerateur / denominateur)
    if denominateur < 0:
        khi += math.pi
    elif numerateur > 0:
        khi += 0.
    else :
        khi += 2 * math.pi
    return khi

def transformation (xin, yin, khi, coef, x0, y0):
    """ Fonction de transformation (xin et yin definit dans [-1,1]) : rotation de khi, puis aggrandissement par coef, puis translation """
    xout = (xin * math.cos (khi) - yin * math.sin (khi)) * coef + x0
    yout = (xin * math.sin (khi) + yin * math.cos (khi)) * coef + y0
    return xout, yout

def object_phase (ax, obs, ra_dec, astre, x0, y0, _diam, couleur_blk, couleur_wht) :
    """ Fonction d'affichage de la phase d'un astre
    ax : un subplot de matplotlib
    obs : objet correspondant a un observateur (lieu et date)
    ra_dec : un booleen True si en (ra,dec) False si en (az,alt)
    astre : objet ephem correspondant a un astre
    x0 et y0 : coordonnees de l'impression dans l'afficheur
    diam : dimension d'affichage de l'astre
    couleur_blk : couleur de la partie dans l'ombre
    couleur_wht : couleur de la partie eclairee
    """
    diam = _diam / 2.
    astre.compute(obs)
    khi = compute_khi(obs, ra_dec, astre)
    if astre.name == 'Moon' :
        k = astre.moon_phase
    else :
        k = astre.phase / 100. # phase en pourcentage
    verts = []
    xpp = []
    ypp = []
    for i in range(-100, 101, 1):
        x = float(i) / 100.
        verts.append([x * diam + x0, math.sqrt(1 - x**2) * diam + y0])      
        if k <= 0.5 :
            yin = math.sqrt((1. - x**2) * (2*k-1)**2)
        else :
            yin = math.sqrt(1. - x**2)
        xout, yout = transformation(x, yin, khi, diam, x0, y0)
        xpp.append(xout)
        ypp.append(yout)
    for i in range(100, -101, -1):
        x = float(i) / 100.
        verts.append([x * diam + x0, -math.sqrt(1 - x**2) * diam + y0])      
        if k > 0.5 :
            yin = -math.sqrt((1. - x**2) * (2*k-1)**2)
        else :
            yin = math.sqrt(1. - x**2)
        xout, yout = transformation(x, yin, khi, diam, x0, y0)
        xpp.append(xout)
        ypp.append(yout)
    verts.append([-diam + x0, y0])
    xout, yout = transformation(-1., 0., khi, diam, x0, y0)
    xpp.append(xout)
    ypp.append(yout)

    polygonBlack = patches.Polygon(verts, color=couleur_blk, lw='1', zorder = '1')
    ax.add_patch(polygonBlack)
    
    verts = []
    for ipt in range(len(xpp)) :
        verts.append([xpp[ipt],ypp[ipt]])      
    polygonWhite = patches.Polygon(verts, color=couleur_wht, lw='1', zorder = '1')
    ax.add_patch(polygonWhite)
    
    #pylab.plot(xpp, ypp, color = couleur_blk, zorder = '1')
    return

def main():
    '''Exemple d'utilisation de la fonction d'affichage d'un astre phase vue depuis la Terre'''
    w,h = pylab.figaspect(1.)
    fig = pylab.figure(figsize = (int(2. * w), int(2. * h)))
    ax = fig.add_subplot(111, axisbg = '0.1') # 0.1 c'est gris fonce
    pylab.xlabel('alt')
    pylab.ylabel('az')
    astre = ephem.Moon()
    obs = ephem.city('Paris')
    obs.date = '2013/3/16 21:000:00'
    astre.compute(obs)
    pylab.title('Moon Phase ('+str(astre.phase)+') - '+str(obs.date))
    # ra_dec = False : Apparent topocentric position represented by the properties : ra and dec, and alt and az.
    object_phase (ax, obs, False, astre, astre.alt * 180. / math.pi, astre.az * 180. / math.pi, astre.size / 1000., '0.2', '0.85')
    pylab.autoscale()
    pylab.show()

if __name__ == "__main__":
    main()
