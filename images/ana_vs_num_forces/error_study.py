#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:13:57 2018

@author: daniel
05/08/2018

Studying error convergence using simple test case. Changing the z position of
x1 and x2, to change the dln line segment and also the proximity to the surface
of the element.

y
^
|
|----->x
 \
  \
   v
    z

x5---------x6
|          |
|          |
|     x1   |
|      \   |
|       \  |
x3-------\-x4
          \
           \
            x2

"""

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator
import scipy.constants as cnst

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
plt.tick_params(axis='y', which='minor')
plt.close('all')
plt.rcParams['lines.markersize'] = 8

def read_parameters(name):
    variables = sio.loadmat("./force_error/%s/%s_params.mat" % (name, name))
    a = variables["a"] = np.asscalar(variables["a"])
    b = variables["b"] = np.ndarray.flatten(variables["b"])
    dist = variables["dist"] = np.ndarray.flatten(variables["dist"])
    mu = variables["mu"] = np.asscalar(variables["mu"])
    n_q = variables["n_q"] = np.ndarray.flatten(variables["n_q"])
    nu = variables["nu"] = np.asscalar(variables["nu"])
    x1i = variables["x1i"] = np.ndarray.flatten(variables["x1i"]) 
    x2i = variables["x2i"] = np.ndarray.flatten(variables["x2i"])
    return a, b, dist, mu, n_q, nu, x1i, x2i

def read_analytic(x1, x2, name, prefix):
    variables = sio.loadmat("./force_error/%s/%s%s_%f_%f.mat" % (name, prefix, name, x1, x2))
    ftota = variables["ftota"] = np.ndarray.flatten(variables["ftota"])
    return ftota

def read_numeric(nq, x1, x2, name, prefix):
    variables = sio.loadmat("./force_error/%s/%s%s_%d_%f_%f.mat" % (name, prefix, name, nq, x1, x2))
    ftotn = variables["ftotn"] = np.ndarray.flatten(variables["ftotn"])
    return ftotn

def read_results(cntr):
    variables = sio.loadmat("./force_error/%s.mat" % str(cntr))
    ftota = variables["ftota"] = np.ndarray.flatten(variables["ftota"])
    ftotn = variables["ftotn"] = np.ndarray.flatten(variables["ftotn"])
    q = variables["q"] = np.ndarray.flatten(variables["q"])
    w = variables["w"] = np.ndarray.flatten(variables["w"])
    x1 = variables["x1"] = np.ndarray.flatten(variables["x1"])
    x2 = variables["x2"] = np.ndarray.flatten(variables["x2"])
    return ftota, ftotn, q, w, x1, x2

def errs(expected, observed, eps = 1E-32):
    diff = observed - expected
    idx = np.abs(diff) < eps
    diff[idx] = 0
    rel_err = np.abs(np.nan_to_num(diff/expected))
    mean_err = np.nan_to_num(np.mean(rel_err))
    rms_rel_err = np.nan_to_num(np.sqrt(np.mean(rel_err**2)))
    return diff, rel_err, mean_err, rms_rel_err

"""
Only forces in y
"""
a, b, dist, mu, nq, nu, x1i, x2i = read_parameters('eperp')

x1 = np.zeros(3)
x2 = np.zeros(3)
a_rel_err = np.ndarray(shape=(1728,3))
a_abs_err = np.ndarray(shape=(1728,3))
dln_len = np.zeros(1728)
idx = np.zeros(145, dtype = int)
distance = np.zeros(144)
quadpts = np.zeros(144, dtype = int)

cntr = -1
cntr2 = 0
for j in np.arange(0, len(dist)):
    x1[2] = x1i[2]
    x1[2] += dist[j]
    for i in np.arange(0, len(nq)):
        cntr2 += 1
        for k in np.arange(0, len(dist)):
            x2[2] = x1[2] + dist[k]
            if x2[2] == x1[2]:
                x2[2] += dist[1]
                continue
            ftota = read_analytic(x1[2], x2[2], 'eperp', 'a')
            ftotn = read_numeric(nq[i], x1[2], x2[2],'eperp', 'n')
            diff, rel_err, mean_err, rms_rel_err = errs(ftota, ftotn, 1E-15)
#            print("%d\t%f\t%f\t%.15f\t%.17f" %(nq[i], x1[2], x2[2]-x1[2],ftotn[1], ftota[1]))
            cntr += 1
            a_rel_err[cntr, :] = rel_err
            a_abs_err[cntr, :] = np.abs(diff)
            dln_len[cntr] = x2[2] - x1[2]
        idx[cntr2] = cntr+1
        distance[cntr2-1] = dist[j]
        quadpts[cntr2-1] = nq[i]
#        print(ftota[1], ftotn[1])
f_lbl = ["x", "y", "z"]
markers=['o', 'v', '8', 's', '^', 'h', '*', 'X', 'D', 'p', 'P', 'H', '<']
cntr = 0
for i in np.arange(0, len(idx)-1):
    mkridx = np.mod(i, len(dist))
    cntr+=1
    if mkridx == 0.0:
        cntr = 0
        f, axarr = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
        f2, axarr2 = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
    if (quadpts[i] > 500 and np.mod(quadpts[i], 10) == 1) or quadpts[i] >= 1500:
        continue
    label = r"%.0f" % distance[i]
    name = r"%f" %distance[i]
    f.suptitle(r"$x^1_z = %s~b$" % label)
    axarr.plot(dln_len[idx[i]: idx[i+1]], a_rel_err[idx[i]:idx[i+1], 1], '-o', label="Q = %d" % quadpts[i], marker = markers[mkridx])
    axarr.set_xscale('log')
    axarr.set_yscale('log')
    axarr.set_ylabel(r"$\log\left|\bm{F}^{\eta}_{%s}\right|$" % f_lbl[1])
    axarr.set_xlabel(r"Dislocation Segment Length, $b$")
    axarr.grid(which="major", linestyle = "-")
    axarr.grid(which="minor", linestyle = "--")
    axarr.legend(loc=1)
#    axarr.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    f.tight_layout()
    f.subplots_adjust(top=0.92)
    f2.suptitle(r"$x^1_z = %s~b$" % label)
    axarr2.plot(dln_len[idx[i]: idx[i+1]], a_abs_err[idx[i]:idx[i+1], 1], '-o', label="Q = %d" % quadpts[i], marker = markers[mkridx])
    axarr2.set_xscale('log')
    axarr2.set_yscale('log')
    axarr2.set_ylabel(r"$\log\left|\bm{F}^{\epsilon}_{%s}\right|$" % f_lbl[1])
    axarr2.set_xlabel(r"Dislocation Segment Length, $b$")
    axarr2.grid(which="major", linestyle = "-")
    axarr2.grid(which="minor", linestyle = "--")
#    axarr2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    f2.tight_layout()
    f2.subplots_adjust(top=0.92)
#    if mkridx == len(markers)-2:
    if cntr == 8:
        f.savefig("perp_e_xz=%s.pdf" % name, format='pdf')
#        f2.savefig("abs_perp_e_xz=%s.pdf" % label, format='pdf')
    
"""
No forces
"""
a, b, dist, mu, nq, nu, x1i, x2i = read_parameters('spar')

x1 = np.zeros(3)
x2 = np.zeros(3)
a_rel_err = np.ndarray(shape=(144,3))
a_abs_err = np.ndarray(shape=(144,3))
dln_len = np.zeros(144)
idx = np.zeros(len(nq)+1, dtype = int)
distance = np.zeros(len(nq)+1)
quadpts = np.zeros(len(nq)+1, dtype = int)

cntr = -1
for i in np.arange(0, len(nq)):
    for j in np.arange(0, len(dist)):
        cntr += 1
        ftota = read_analytic(dist[j], dist[j], 'spar', 'a')
        ftotn = read_numeric(nq[i], dist[j], dist[j],'spar', 'n')
#        print(dist[j], ftota, ftotn)
        diff, rel_err, mean_err, rms_rel_err = errs(ftota, ftotn, 1E-15)
        a_rel_err[cntr, :] = rel_err
        a_abs_err[cntr, :] = np.abs(diff)
        dln_len[cntr] = dist[j]
    idx[i+1] = cntr+1
    quadpts[i] = nq[i]

f_lbl = ["x", "y", "z"]
markers=['o', 'v', '8', 's', '^', 'h', '*', 'X', 'D', 'p', 'P', 'H', '<']
for j in np.arange(0,3):
    for i in np.arange(0, len(idx)-1):
        mkridx = np.mod(i, len(dist))
        if mkridx == 0.0:
            f, axarr = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
            f2, axarr2 = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
        axarr.plot(dln_len[idx[i]: idx[i+1]], a_rel_err[idx[i]:idx[i+1], j], '-o', label="Q = %d" % quadpts[i], marker = markers[mkridx])
        axarr.set_xscale('log')
        axarr.set_yscale('log')
        axarr.set_ylabel(r"$\log\left|\bm{F}^{\eta}_{%s}\right|$" % f_lbl[j])
        axarr.set_xlabel(r"Dislocation Distance from SE, $b$")
        axarr.grid(which="major", linestyle = "-")
        axarr.grid(which="minor", linestyle = "--")
#        axarr.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        f.tight_layout()
        f.subplots_adjust(top=0.92)
        axarr2.plot(dln_len[idx[i]: idx[i+1]], a_abs_err[idx[i]:idx[i+1], j], '-o', label="Q = %d" % quadpts[i], marker = markers[mkridx])
        axarr2.set_xscale('log')
        axarr2.set_yscale('log')
        axarr2.set_ylabel(r"$\log\left|\bm{F}^{\epsilon}_{%s}\right|$" % f_lbl[j])
        axarr2.set_xlabel(r"Dislocation Distance from SE, $b$")
        axarr2.grid(which="major", linestyle = "-")
        axarr2.grid(which="minor", linestyle = "--")
#        axarr2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        f2.tight_layout()
        f2.subplots_adjust(top=0.92)
    #    if mkridx == len(markers)-2:
    #        f.savefig("par_s.pdf", format='pdf')
    #        f2.savefig("abs_par_s.pdf", format='pdf')
        
"""
The z coordinates have forces.
"""
a, b, dist, mu, nq, nu, x1i, x2i = read_parameters('epar')

x1 = np.zeros(3)
x2 = np.zeros(3)
a_rel_err = np.ndarray(shape=(144,3))
a_abs_err = np.ndarray(shape=(144,3))
dln_len = np.zeros(144)
idx = np.zeros(len(nq)+1, dtype = int)
distance = np.zeros(len(nq)+1)
quadpts = np.zeros(len(nq)+1, dtype = int)

cntr = -1
for i in np.arange(0, len(nq)):
    for j in np.arange(0, len(dist)):
        cntr += 1
        ftota = read_analytic(dist[j], dist[j], 'epar', 'a')
        ftotn = read_numeric(nq[i], dist[j], dist[j],'epar', 'n')
#        print(dist[j], ftota, ftotn)
        diff, rel_err, mean_err, rms_rel_err = errs(ftota, ftotn, 1E-15)
        a_rel_err[cntr, :] = rel_err
        a_abs_err[cntr, :] = np.abs(diff)
#        if np.abs(diff[1]) > 10:
#            print(ftota,ftotn)
        dln_len[cntr] = dist[j]
    idx[i+1] = cntr+1
    quadpts[i] = nq[i]
    
f_lbl = ["x", "y", "z"]
markers=['o', 'v', '8', 's', '^', 'h', '*', 'X', 'D', 'p', 'P', 'H', '<']
cntr = 0
for j in np.arange(2,3):
    for i in np.arange(0, len(idx)-1):
        mkridx = np.mod(i, len(dist))
        cntr +=1
        if mkridx == 0.0:
            cntr = 0
            f, axarr = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
            f2, axarr2 = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
        if (quadpts[i] > 500 and np.mod(quadpts[i], 10) == 1) or quadpts[i] >= 1500:
            continue
        axarr.plot(dln_len[idx[i]: idx[i+1]], a_rel_err[idx[i]:idx[i+1], j], '-o', label="Q = %d" % quadpts[i], marker = markers[mkridx])
        axarr.set_xscale('log')
        axarr.set_yscale('log')
        axarr.set_ylabel(r"$\log\left|\bm{F}^{\eta}_{%s}\right|$" % f_lbl[j])
        axarr.set_xlabel(r"Dislocation Distance from SE, $b$")
        axarr.grid(which="major", linestyle = "-")
        axarr.grid(which="minor", linestyle = "--")
        axarr.legend(loc=1)
#        axarr.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        f.tight_layout()
        f.subplots_adjust(top=0.92)
        axarr2.plot(dln_len[idx[i]: idx[i+1]], a_abs_err[idx[i]:idx[i+1], j], '-o', label="Q = %d" % quadpts[i], marker = markers[mkridx])
        axarr2.set_xscale('log')
        axarr2.set_yscale('log')
        axarr2.set_ylabel(r"$\log\left|\bm{F}^{\epsilon}_{%s}\right|$" % f_lbl[j])
        axarr2.set_xlabel(r"Dislocation Distance from SE, $b$")
        axarr2.grid(which="major", linestyle = "-")
        axarr2.grid(which="minor", linestyle = "--")
#        axarr2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        f2.tight_layout()
        f2.subplots_adjust(top=0.92)
#        if mkridx == len(markers)-2:
        if cntr == 8:
            f.savefig("par_e_%s.pdf"%f_lbl[j], format='pdf')
#            f2.savefig("abs_par_e_%s.pdf"%f_lbl[j], format='pdf')






#for i in np.arange(0, len(idx)-1):
#    f, axarr = plt.subplots(1, sharex = True, figsize=(7,7/cnst.golden))
#    if distance[i] == 0:
#        label = "0.0"
#    else:
#        label = r"10^{%.0f}/2" % np.log10(distance[i])
#    f.suptitle(r"$x^1_z = %s$" % label)
#    for j in np.arange(0,1):
#        # if more than one subplot use axarr[j], and dln_len[foo:bar,wig:hat,j]
#        axarr.plot(dln_len[idx[i]: idx[i+1]], a_rel_err[idx[i]:idx[i+1], j+1], '-o', label="Q = %d" % quadpts[i])
#    plt.xscale('log')
#    plt.ylabel(r"$\bm{F}^{\eta}_{%s}$" % f_lbl[j+1])
##        if j == 2:
##            axarr[j].set_xlabel(r"Dislocation Segment Length")
#    plt.xlabel(r"Dislocation Segment Length, $2/\mathrm{d}x$")
#    plt.grid(which="major", linestyle = "-")
#    plt.grid(which="minor", linestyle = "--")
#    plt.legend(loc=0)
#    plt.tight_layout()
#    f.subplots_adjust(top=0.92)
#    if quadpts[i] == 10 and np.log10(distance[i]) == -1:
#        plt.savefig("xz=10_q=10.pdf")
#    if quadpts[i] == 11 and np.log10(distance[i]) == -1:
#        plt.savefig("xz=10_q=11.pdf")
##    f.savefig("xz=%s.pdf" % label)
            
