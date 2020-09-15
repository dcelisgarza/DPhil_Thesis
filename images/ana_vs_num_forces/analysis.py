# -*- coding: utf-8 -*-
"""
Daniel Celis Garza
18/07/2018

Analysis of .mat Files of numerical vs analytical surface forces.
"""

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 25)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
plt.close('all')

def read_mat(mx, face):
    variables = sio.loadmat("mx=%s_face=%s.mat" % (mx, face))
    dim = np.ndarray.flatten(variables["dim"])
    f_dln_a_lin = variables["ftilda_a_lin"]
    f_dln_a_par = variables["ftilda_a_par"]
    f_dln_n = variables["ftilda_n"]
    xmid = variables["midpoint_element"]
    return dim, f_dln_a_lin, f_dln_a_par, f_dln_n, xmid

def select_face(face, dim, xmid):
    # min(y) xz plane
    if face == 1:
        X = np.reshape(xmid   [:, 0], (dim[2], dim[1]))
        Y = np.reshape(xmid   [:, 2], (dim[2], dim[1]))
        lbl = ['x', 'z']
    # max(z) xy plane
    elif face == 2:
        X = np.reshape(xmid   [:, 0], (dim[2], dim[1]))
        Y = np.reshape(xmid   [:, 1], (dim[2], dim[1]))
        lbl = ['x', 'y']
    # max(x) yz plane
    elif face == 3:
        X = np.reshape(xmid   [:, 1], (dim[2], dim[1]))
        Y = np.reshape(xmid   [:, 2], (dim[2], dim[1]))
        lbl = ['y', 'z']
    return X, Y, lbl

def errs(expected, observed, eps = 1E-13):
    diff = observed - expected
    idx = np.abs(diff) < eps
    diff[idx] = 0
    rel_err = np.nan_to_num(np.abs(np.abs(diff)/expected))
    mean_err = np.array(np.mean(rel_err, 0))
    rms_rel_err = np.array(np.sqrt(np.mean(rel_err**2, 0)))
    diff = np.abs(diff)
    return rel_err, mean_err, rms_rel_err, diff

def plot_figs(dim, lbl, X, Y, Zi, z_lbl, figname, levels = 10, cmap = cm.terrain, norm_arr = np.zeros((2,3))):
    f, axarr = plt.subplots(3)
    f.set_figheight(10)
    f.set_figwidth(10)
    f_lbl = ['x', 'y', 'z']
    for i in np.arange(0,3,1):
        Z = np.reshape(Zi[:, i], (dim[2], dim[1]))
        if np.all(norm_arr == 0):
            norm = cm.colors.Normalize(vmax=np.max(Z), vmin=np.min(Z))
        else:
            norm = cm.colors.Normalize(vmax=np.max(norm_arr[:,i]), vmin=np.min(norm_arr[:,i]))
#        filled = axarr[i].contourf(X,Y,Z,levels,cmap=cm.get_cmap(cmap, levels),norm=norm)
        if dim[2] == dim[1]:
            Z = Z.transpose()
        filled = axarr[i].imshow(Z,interpolation='bicubic',extent=[X.min(), X.max(), Y.min(), Y.max()],cmap=cmap,norm=norm, origin='lower')
        cbar = plt.colorbar(filled, ax = axarr[i], format = '%1.3e')
#        cbar.set_label(r'$\sigma_{%s}\,\mu^{-1}$'%lbl[0], rotation = 90)
#        if dim[2] == dim[1]:
#            Z = Z.transpose()
#        lines = axarr[i].contour(X,Y,Z,levels,linestyles='-')
#        plt.colorbar(lines, ax = axarr[i])
        axarr[i].set_title (r"$\bm{F}^{%s}_{%s}$" % (z_lbl, f_lbl[i]))
        axarr[i].set_xlabel(r"$%s,~b$" % lbl[0])
        axarr[i].set_ylabel(r"$%s,~b$" % lbl[1])
        axarr[i].set_aspect('equal')
        axarr[i].set_xticks(np.arange(0,50000,10000))

#        axarr[i].set_xticks(np.arange(0, 50000, 10000));
#        axarr[i].set_yticks(np.arange(0, 15000, 5000));
    plt.tight_layout()
    
    plt.savefig('%s_%s' % (z_lbl, figname), bbox_inches='tight')
    return

mx_min = 60#10
mx_max = 120#120
mx_stp = 60##10
face_min = 1
face_max = 1
face_stp = 1
eps = 1E-9

face_arr = np.arange(face_min, face_max + face_stp, face_stp)
mx_arr = np.arange(mx_min, mx_max + mx_stp, mx_stp)
levels = 60

cmap = cm.terrain

n_plots = int((mx_max - mx_min + mx_stp)/mx_stp)
rms_rel_err_arr = np.zeros([face_max, n_plots, 3]);

for face in face_arr:
    for mx in mx_arr:
        dim, f_dln_a_lin, f_dln_a_par, f_dln_n, xmid = read_mat(mx, face)
        X, Y, lbl = select_face(face, dim, xmid)
        
        diff = f_dln_a_lin-f_dln_a_par
        
        # Parallel
#        rel_error, mean_err, rms_rel_err, diff = errs(f_dln_a_par, f_dln_n, eps)
#        rms_rel_err_arr[int((face-face_min)/face_stp - 1), int((mx-mx_min)/mx_stp - 1), :] = rms_rel_err
#        norm_arr = np.array([np.min(f_dln_n, 0)/(mean_err+1),
#                             np.max(f_dln_n, 0)/(mean_err+1)])
#        plot_figs(dim, lbl, X, Y, f_dln_a_par, 'a_p', 'mx=%s_face=%s.pdf' % (mx, face), levels = levels, 
#                  cmap = cmap, norm_arr = norm_arr)
#        
#        # Error parallel vs serial
#        rel_error, mean_err, rms_rel_err, diff = errs(f_dln_a_par, f_dln_a_lin, eps)
#        rms_rel_err_arr[int((face-face_min)/face_stp - 1), int((mx-mx_min)/mx_stp - 1), :] = rms_rel_err
#        norm_arr = np.array([np.min(rel_error, 0), np.max(rel_error, 0)])
#        plot_figs(dim, lbl, X, Y, rel_error, '\eta_{lp}', 'mx=%s_face=%s.pdf' % (mx, face), levels = levels, 
#                  cmap = cmap, norm_arr = norm_arr)
        
        # Serial
        rel_error, mean_err, rms_rel_err, diff = errs(f_dln_a_lin, f_dln_n, eps)
#        rms_rel_err_arr[int((face-face_min)/face_stp - 1), int((mx-mx_min)/mx_stp - 1), :] = rms_rel_err
#        norm_arr = np.array([np.min(f_dln_n, 0)/(mean_err+1),
#                             np.max(f_dln_n, 0)/(mean_err+1)])
#        plot_figs(dim, lbl, X, Y, f_dln_a_lin, 'a_l', 'mx=%s_face=%s.pdf' % (mx, face), levels = levels, 
#                  cmap = cmap, norm_arr = norm_arr)
#        
#        # Numeric
#        norm_arr = np.array([np.min(f_dln_a_lin, 0)*(mean_err+1),
#                             np.max(f_dln_a_lin, 0)*(mean_err+1)])
#        plot_figs(dim, lbl, X, Y, f_dln_n, 'n', 'mx=%s_face=%s.pdf' % (mx, face), levels = levels,
#                  cmap = cmap, norm_arr = norm_arr)
        
        # Error numeric vs analytic
        if mx == mx_min:
            norm_arr = np.array([np.min(rel_error, 0), np.max(rel_error, 0)])
        else:
            norm_arr = norm_arr
#        norm_arr = np.array([np.min(diff, 0), np.max(diff, 0)])
        plot_figs(dim, lbl, X, Y, rel_error, '\eta', 'mx=%s_face=%s.pdf' % (mx, face), levels = levels,
                  cmap = cmap, norm_arr = norm_arr)
#        print(norm_arr)

#x = np.arange(mx_min, mx_max + mx_stp, mx_stp)
#
#for face in face_arr:
#    y1 = rms_rel_err_arr[face-1, : , 0]
#    y2 = rms_rel_err_arr[face-1, : , 1]
#    y3 = rms_rel_err_arr[face-1, : , 2]
#    fig = plt.figure()
#    plt.plot(x, y1, label = r'$F_x$')
#    plt.plot(x, y2, label = r'$F_y$')
#    plt.plot(x, y3, label = r'$F_z$')
#    plt.xlabel(r'Nodes in x-dimension')
#    plt.ylabel(r'$\mathrm{RMS}(\eta)$')
#    plt.legend(loc=0)
#    plt.tight_layout()
#    plt.savefig('rms_%s.pdf' % face, bbox_inches='tight')
    