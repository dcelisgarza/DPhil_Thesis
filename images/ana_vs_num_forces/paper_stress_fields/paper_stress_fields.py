#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:21:11 2018

@author: daniel
"""
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 25)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
plt.close('all')

def read_parameters(dln_type, dln_orientation):
    variables = sio.loadmat("%s%s_params.mat" % (dln_type, dln_orientation))
    a = variables["a"] = np.asscalar(variables["a"])
    b = variables["b"] = np.ndarray.flatten(variables["b"])
    dist = variables["dist"] = np.ndarray.flatten(variables["dist"])
    x1 = variables["x1"] = np.ndarray.flatten(variables["x1"])
    x2 = variables["x2"] = np.ndarray.flatten(variables["x2"])
    mu = variables["mu"] = np.asscalar(variables["mu"])
    nu = variables["nu"] = np.asscalar(variables["nu"])
    return a, b, dist, x1, x2, nu, mu

def read_fields(dln_type, dln_orientation, dist):
    variables = sio.loadmat("%s%s_%f.mat" % (dln_type, dln_orientation, dist))
    sxx = variables['SXX']
    syy = variables['SYY']
    szz = variables['SZZ']
    sxy = variables['SXY']
    sxz = variables['SXZ']
    syz = variables['SYZ']
    X = variables['X']
    Y = variables['Y']
    return [X, Y], [sxx, syy, szz, sxy, sxz, syz]

def plot_figs(lbl, X, Y, Z, levels = 20, cmap = cm.terrain, norm = np.zeros((2,1)), figname = ''):
    f, ax = plt.subplots(1)
    f.set_figheight(10)
    f.set_figwidth(10)
    if norm.all() == 0:
        norm = cm.colors.Normalize(vmax=np.max(Z), vmin=np.min(Z))
    else:
        norm = cm.colors.Normalize(vmax=np.max(norm), vmin=np.min(norm))
    filled = ax.contourf(X,Y,Z,levels,cmap=cm.get_cmap(cmap, levels),norm=norm)
#    filled = ax.imshow(Z,interpolation='bilinear',extent=[X.min(), X.max(), Y.min(), Y.max()],cmap=cmap,norm=norm, origin='lower')
    cbar = plt.colorbar(filled, ax = ax, format = r'%1.0e',fraction=0.046, pad=0.04)
    cbar.set_label(r'$\sigma_{%s}\,\mu^{-1}$'%lbl[0], rotation = 90)
    ax.set_title (r"$%s,~R = %.0f\, b$" % (lbl[1], lbl[2]))
    ax.set_xlabel(r"$x,\,b$")
    ax.set_ylabel(r"$y,\,b$")
    ax.set_aspect('equal')
    plt.tight_layout()
    if figname != '':
        plt.savefig('%s_%s%s_%f%s' % (lbl[0], lbl[3], lbl[4], lbl[2], figname), bbox_inches='tight')
    return

dln_type = 'e'
dln_orientation = 'par'
stress_labels = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']

a, b, dist, x1, x2, mu, nu = read_parameters(dln_type, dln_orientation)

if dln_orientation == 'perp':
    auxorien = 'perp'
elif dln_orientation == 'par':
    auxorien = 'parallel'
if dln_type == 's':
    auxtype = '\mathrm{screw}'
elif dln_type == 'e':
    auxtype = '\mathrm{edge}'

#dist = np.ones(1)
#dist = 5*dist
#dln_type = 't'+dln_type
dln_orientation = 'par_pm4'
for i in np.arange(len(dist)):
    r = dist[i]
    XY, stresses = read_fields(dln_type, dln_orientation, r)
    for j in np.arange(len(stress_labels)):
        lbl = ['%s' % stress_labels[j], '\%s{}~%s' % (auxorien, auxtype), r, dln_type, dln_orientation]
        try:
            plot_figs(lbl, XY[0], XY[1], stresses[j], figname='.pdf')
        except:
            pass
