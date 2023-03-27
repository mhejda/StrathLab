# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:20:11 2023

@author: Dafydd Owen-Newns
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def plot_square(x, y):
    x_plot = []
    y_plot = []
    for i in range(len(x)-1):
        x_plot.append(x[i])
        x_plot.append(x[i+1])
        y_plot.append(y[i])
        y_plot.append(y[i])
    x_plot.append(x[-1])
    y_plot.append(y[-1])
    fig, ax = plt.subplots()
    ax.plot(x_plot, y_plot)
    
def square_wave(time_up, time_down, sr = 12e9):

    up_len = int(sr*time_up)
    down_len = int(sr*time_down)
    
    wf = np.concatenate((np.ones(up_len), np.zeros(down_len)))
    
    return wf
    
def create_marker(tail = 10e-9, sr = 10e9):
    s = square_wave(0.5e-9, 0.5e-9, sr)
    return np.concatenate((s, s, s, np.zeros(int(tail*sr))))

def mesh(M, T = None, size = (8, 4.8)):
    fig, ax = plt.subplots(1, 1, figsize = size)
    if(T is None):
        ax.pcolormesh(M)
    else:
        ax.pcolormesh(M, cmap = T)
    return fig, ax

def draw(x, y, size=(8, 4.8)):
    fig, ax = plt.subplots(1, 1, figsize = size)
    ax.plot(x, y)
    return fig, ax

def draw_y(y, size=(8, 4.8)):
    fig, ax = plt.subplots(1, 1, figsize = size)
    ax.plot(range(len(y)), y)
    return fig, ax