#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 13:35:00 2018

@author: daniel
"""

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
plt.close('all')

def read_parameters(dln_type, dln_orientation):
    variables = sio.loadmat("./stress_fields/%s%s_field/%s%s_params.mat" % (dln_type, dln_orientation, dln_type, dln_orientation))
    a = variables["a"] = np.asscalar(variables["a"])
    b = variables["b"] = np.ndarray.flatten(variables["b"])
    dist = variables["dist"] = np.ndarray.flatten(variables["dist"])
    x1 = variables["x1"] = np.ndarray.flatten(variables["x1"])
    x2 = variables["x2"] = np.ndarray.flatten(variables["x2"])
    mu = variables["mu"] = np.asscalar(variables["mu"])
    nu = variables["nu"] = np.asscalar(variables["nu"])
    return a, b, dist, x1, x2, nu, mu

def read_fields(dln_type, dln_orientation, dist):
    variables = sio.loadmat("./stress_fields/%s%s_field/%s%s_%f.mat" % (dln_type, dln_orientation, dln_type, dln_orientation, dist))
    sxx = variables['SXX']
    syy = variables['SYY']
    szz = variables['SZZ']
    sxy = variables['SXY']
    sxz = variables['SXZ']
    syz = variables['SYZ']
    return [sxx, syy, szz, sxy, sxz, syz]

def plot_figs(lbl, X, Y, Z, levels = 50, cmap = cm.terrain, norm = np.zeros((2,1)), figname = ''):
    f, ax = plt.subplots(1)
    f.set_figheight(10)
    f.set_figwidth(10)
    if norm.all() == 0:
        norm = cm.colors.Normalize(vmax=np.max(Z), vmin=np.min(Z))
    else:
        norm = cm.colors.Normalize(vmax=np.max(norm), vmin=np.min(norm))
    filled = ax.contourf(X,Y,Z,levels,cmap=cm.get_cmap(cmap, levels),norm=norm)
#    filled = ax.imshow(Z,interpolation='bicubic',extent=[X.min(), X.max(), Y.min(), Y.max()],cmap=cmap,norm=norm, origin='lower')
    plt.colorbar(filled, ax = ax, format = '%2.3f')
    ax.set_title (r"$\bm{\sigma_{%s}},~%s,~R = %f$" % (lbl[0], lbl[1], lbl[2]))
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_aspect('equal')
    plt.tight_layout()
    if figname != '':
        plt.savefig('%s_%s%s_%f%s' % (lbl[0], lbl[3], lbl[4], lbl[2], figname), bbox_inches='tight')
    return


x = np.linspace(-1,1,101)
y = x
X, Y = np.meshgrid(x,y)

r = 0.
dln_type = 'e'
dln_orientation = 'perp'
stress_labels = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']

a, b, dist, x1, x2, nu, mu = read_parameters(dln_type, dln_orientation)

if dln_orientation == 'perp':
    auxorien = 'perp'
elif dln_orientation == 'par':
    auxorien = 'parallel'
if dln_type == 's':
    auxtype = '\mathrm{screw}'
elif dln_type == 'e':
    auxtype = '\mathrm{edge}'

for i in np.arange(len(dist)):
    r = dist[i]
    stresses = read_fields(dln_type, dln_orientation, r)
    for j in np.arange(len(stress_labels)):
        lbl = ['%s' % stress_labels[j], '\%s{}~%s' % (auxorien, auxtype), r, dln_type, dln_orientation]
        try:
            plot_figs(lbl, X, Y, stresses[j], figname='.pdf')
        except:
            pass

