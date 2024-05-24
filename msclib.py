# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:50:35 2023

@author: David Fox
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import ipywidgets as widgets
from scipy.optimize import fsolve
from scipy import linalg
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt



# Find roots of a function from R^2 to R^2
def find_roots(f, xrange, yrange, n, args, tol=1e-8):
    x1, x2 = xrange
    y1, y2 = yrange
    x = np.linspace(x1, x2, n)
    y = np.linspace(y1, y2, n)
    xv, yv = np.meshgrid(x, y)
    
    roots = []
    
    for i in range(n):
        for j in range(n):
            root = fsolve(f, [xv[i][j], yv[i][j]], args, maxfev=9000)
            if np.linalg.norm(f(root, *args)) < tol:
                if len(roots) == 0:
                    roots.append(root)
                elif all([not(np.allclose(x, root)) for x in roots]):
                    roots.append(root)
                
    return roots

# Plotting Phase portraits
def color_map(x):
    colors = np.log10(x)
    return np.nan_to_num(colors, neginf=0) 
    
def phase_portrait(ax, f, t, xlim, ylim, num_pts, args, norm_arrows=True, stream=True, quiver=True):
    
    x = np.linspace(-xlim, xlim, num_pts)
    y = np.linspace(-ylim, ylim, num_pts)
    X, Y = np.meshgrid(x, y)
    u = np.zeros((num_pts, num_pts))
    v = np.zeros((num_pts, num_pts))
    flow_mag = np.zeros((num_pts, num_pts))

    
    for i in range(num_pts):
        for j in range(num_pts):
            u[i,j], v[i,j] = f(0, [X[i,j], Y[i,j]], *args)
            flow_mag[i,j] = np.sqrt(u[i,j]**2 + v[i,j]**2)
            if (u[i,j], v[i,j]) != (0,0):
                u[i,j] *= 1/flow_mag[i,j]
                v[i,j] *= 1/flow_mag[i,j]
    
    if quiver:
        ax.quiver(X[::1, ::1], Y[::1, ::1], u[::1, ::1], v[::1, ::1], color_map(flow_mag[::1, ::1]), cmap="winter", alpha=0.1)
        
    ax.streamplot(X, Y, u, v,  color=color_map(flow_mag), cmap="winter", density=1)
    ax.grid(True, which='both') 


# Classify a fixed point of a 2d dynamical system
def classify_fp(jac, x, params):
    evals = linalg.eigvals(jac(x, *params))
    
    # Check if fixed point is a stable or unstable.
    if all(evals.real < 0):
        stability = 'stable'
    else:
        stability = 'unstable'
    
    # Check if fixed point is a spiral, node or saddle.
    if not(all(evals.imag == 0)):
        character = 'spiral'
    elif all(evals < 0) or all(evals > 0):
        character = 'node'
    else:
        character = 'saddle'
    
    return (stability, character)


# Equation used to plot streamlines
def f_2(t, r, gamma, m11, m12, m21, m22, u1, u2):
    r1, r2 = r
    r_next = np.zeros(2)
    r_next[0] = gamma*(-r1 + np.tanh(m11*r1 + m12*r2 + u1))
    r_next[1] = gamma*(-r2 + np.tanh(m21*r1 + m22*r2 + u2))
    return r_next

# Equations for r1, r2 nullclines in 2D system with fixed mij. Used for plotting nullclines
def r1_null_2(r1, r2, m11, m12, u1):
    return -r1 + np.tanh(m11*r1 + m12*r2 + u1)

def r2_null_2(r1, r2, m21, m22, u2):
    return -r2 + np.tanh(m21*r1 + m22*r2 + u2)

# Equations for nullclines used to find intersections (fixed points)
def nullclines_2(r, m11, m12, m21, m22, u1, u2):
    r1, r2 = r
    return np.asarray([r1_null_2(r1, r2, m11, m12, u1), r2_null_2(r1, r2, m21, m22, u2)])

def jac(r, gamma, m11, m12, m21, m22, u1, u2):
    r1, r2 = r
    j11 = gamma*(-1 + m11*1/np.cosh(m11*r1 + m12*r2 + u1)**2)
    j12 = gamma*m12*1/np.cosh(m11*r1 + m12*r2 + u1)**2
    j21 = gamma*m21*1/np.cosh(m21*r1 + m22*r2 + u2)**2
    j22 = gamma*(-1 + m22*1/np.cosh(m21*r1 + m22*r2 + u2)**2)
    return np.asarray([[j11, j12], [j21, j22]])


def CTRNN_phase_portrait(gamma=0.1, m11=0.1, m12=0.1, m21=0.1, m22=0.1, u1=0.1, u2=0.1, ax=None):    
    if ax == None:
        fig, ax = plt.subplots(1, 1, figsize=(7,7))
    
    # Create phase portrait for (r1, r2) subsystem with fixed values of weights.
    
    # Plot r1, r2 nullclines for fixed m values and (r1, r2) trajectory
    delta = 0.025
    x = np.arange(-2, 2, delta)
    y = np.arange(-2, 2, delta)
    p, q = np.meshgrid(x, y)
    z1 = r1_null_2(p, q, m11, m12, u1)
    z2 = r2_null_2(p, q, m21, m22, u2)
    ax.contour(p, q, z1, [0], colors=["y"])
    ax.contour(p, q, z2, [0], colors=["r"])
    ax.plot([], [], 'r', label='$\dot{r}_1 = 0$')
    ax.plot([], [], 'y', label='$\dot{r}_2 = 0$')
    
    # Find Fixed points of 2D (r1, r2) system by finding intersection of r nullclines.
    params = (gamma, m11, m12, m21, m22, u1, u2)
    fixed_points = find_roots(nullclines_2, [-1, 1], [-1, 1], 20, (params[1:]))
    
    for x in fixed_points:
        fp_type = classify_fp(jac, x, params)
        marker = {'node':'ks', 'spiral':'ko', 'saddle':'k^'}[fp_type[1]]
        if fp_type[0] == 'stable':
            ax.plot(x[0], x[1], marker, markersize=10)
        else:
            ax.plot(x[0], x[1], marker, markerfacecolor='None', markersize=10)
        
    
    # Plot phase portrait for 2D (r1, r2) system
    r_range = np.linspace(-0.99, 0.99, 10000)
    phase_portrait(ax, f_2, np.linspace(0, 500, 20000), 2, 2, 15, params, quiver=False)
    
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel('$r_{1}$', fontsize=20)
    ax.set_ylabel('$r_{2}$', fontsize=20)
    
    return ax
    

def slow_manifold(m_range, m11, m22, u1, u2, alpha):
    # Compute all fixed points of the fast subsystem which make up the slow manifold.
    points = []
    for m in m_range:
        roots = find_roots(nullclines_2, [-1, 1], [-1, 1], 20, (m11, m, alpha*m, m22, u1, u2))
        [points.append(np.append(root, m)) for root in roots]
    
    stable_points = []
    unstable_points = []

    # Compute stability of the points w.r.t the fast subsystem dynamics, then
    # split points into stable and unstable.
    for x in points:
        params = (1, m11, x[2], alpha*x[2], m22, u1, u2)
        fp_type = classify_fp(jac, x[0:2], params)
        if fp_type[0] == 'stable':
            stable_points.append(x)
        elif fp_type[0] == 'unstable':
            unstable_points.append(x)
            
    return (stable_points, unstable_points)


# I bifurcation diagram for 1 neuron as a function of r, paramaterised by self weight m.
def I(r, m):
    return 1/2*np.log((1+r)/(1-r)) - m*r


def f(r, m, I):
    return -r + np.tanh(m*r + I)

def f_prime(r, m, I):
    return -1 + m*(1 - np.tanh(m*r + I)**2)

def I(r, m):
    return 1/2*np.log((1+r)/(1-r)) - m*r

def r_plus(m):
    return np.sqrt((m-1)/m)

def I_bif_curve(ax, r, m):
    if m < 1:
        ax.plot(I(r, m), r, 'g', linewidth=4)
    else:
        left_branch = r[r<-r_plus(m)]
        middle_branch = r[np.logical_and(r>=-r_plus(m),  r<=r_plus(m))]
        right_branch = r[r>r_plus(m)]
        ax.plot(I(left_branch, m), left_branch, 'g', linewidth=4)
        ax.plot(I(middle_branch, m), middle_branch, 'r--', linewidth=4)
        ax.plot(I(right_branch, m), right_branch, 'g', linewidth=4)
        
def I_bif_diag(m):
    r = np.linspace(-1+1e-6, 1-1e-6, 10000)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    I_bif_curve(ax, r, m)
    
    if m > 1:
        r_p = r_plus(m)
        ax.plot(I(r_p, m), r_p, 'ko', markersize=12)
        ax.plot(I(-r_p, m), -r_p, 'ko', markersize=12)
        ax.text(I(r_p, m)-0.7, r_p-0.03, r'$\left(I_{-}, r_{-}\right)$', fontsize=30)
        ax.text(I(-r_p, m)+0.08, -r_p - 0.03, r'$\left(I_{+}, r_{+}\right)$', fontsize=30)
    
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    
    ax.tick_params(direction='out', length=6, width=2)

    xticks = [-1, 0, 1]
    plt.xticks(fontsize=30)
    ax.set_xticks(xticks)
    plt.yticks(fontsize=30)
    ax.set_yticks(xticks)
    
    ax.set_xlabel(r'$I$', fontsize=30)
    ax.set_ylabel(r'$r^{*}$', fontsize=30, rotation=0)
    ax.set_aspect('equal', adjustable='box')
    
    return (fig, ax)
        