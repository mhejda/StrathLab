# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:12:13 2021

@author: thisi
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from tkinter import *
from tkinter import filedialog
from math import sqrt, ceil

plt.rcParams['font.family'] = 'sans-serif'
plt.style.use('default')

def factor_int_2(n):
    #source: https://stackoverflow.com/questions/39248245/factor-an-integer-to-something-as-close-to-a-square-as-possible
    val = ceil(sqrt(n))
    while True:
        if not n%val:
            val2 = n//val
            break
        val -= 1
    return val, val2

def MovingAverage(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def PlotGradient(ax_object,curvex, curvey,color='xkcd:magenta',ymin=0,ymax=0.02):
    # Plot spike with gradient
    line, = ax_object.plot(curvex, curvey, color=color, alpha=1.0,lw=1)

    fill_color = line.get_color()
    
    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 0.75 if alpha is None else alpha
    
    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:,:,:3] = rgb
    z[:,:,-1] = np.linspace(0, alpha, 100)[:,None]
    
    xmin, xmax = curvex.min(), curvex.max()  #curvey.min, curvey.max()
    im = ax_object.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)
    
    xy = np.column_stack([curvex, curvey])
    xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    clip_path = Polygon(xy, facecolor='none', edgecolor='none', closed=True)
    ax_object.add_patch(clip_path)
    im.set_clip_path(clip_path)  

def PlotSingle(fname, fig, ax,alpha = 1.0,labelled_x = True):
    axis_upscale=1e9  
    MA_len = 50
    loaded = np.load(fname)
    xx = loaded[:,0]
    yy = loaded[:,1]
    dt = (xx[100]-xx[0])/100

    yyMA = MovingAverage(yy,MA_len)
    xxMA = xx[0:yyMA.shape[0]]                
    xxFrom0 = np.arange(0,dt*xxMA.shape[0],dt)

    ax.plot(xxFrom0*axis_upscale, yyMA,color='xkcd:crimson',label=fname[0:15],alpha=alpha)
    if labelled_x:
        ax.set_xlabel(r't [ns]')
    #ax.set_ylabel(r'OSC readout [V]')
    ax.legend(loc=1,fontsize='x-small')

def QuickBatchRenamer(set_folder):    
    os.chdir(set_folder)
    _, _, filenames = next(os.walk(os.getcwd()), (None, None, []))
    for file in filenames:
        if file[-4:] == '.npy':
            newfn = file.replace("SubSamp", "ReSamp")
            os.rename(file, newfn)        

def CrunchCSVs(set_folder,del_csvs=True):
    os.chdir(set_folder)
    
    source = 10
    subsample = 1
    appendix = " "
    
    _, _, filenames = next(os.walk(os.getcwd()), (None, None, []))
    for file in filenames:
        if file[-4:] == '.csv' and '{Y}' in file:
            
            
            if '2.5TSa' in file:
                source=2500
                subsample = 250
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 2.5TSa    ",end='')
            elif '2TSa' in file:
                source=2000
                subsample = 200 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 2 TSa   ",end='')    
            elif '1TSa' in file:
                source=1000
                subsample = 100 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 1 TSa   ",end='')    
            elif '500GSa' in file:
                source=250
                subsample = 25 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 500 GSa   ",end='')                    
            elif '500GSa' in file:
                source=250
                subsample = 25 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 500 GSa   ",end='')    
                
            elif '250GSa' in file:
                source=250
                subsample = 25 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 250 GSa   ",end='')
            elif '100GSa' in file:
                source=100
                subsample = 10 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 100 GSa   ",end='')   
                
            elif '50GSa' in file:
                source=50
                subsample = 5 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 50 GSa   ",end='')                   
            
            elif '10GSa' in file:
                source=10
                subsample = 1 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 10 GSa   ",end='')        
            elif '1GSa' in file:
                source=1
                subsample = 1 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 1 GSa   ",end='')
                
            elif '500MSa' in file:
                source=0.5
                subsample = 1 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 500 MSa   ",end='')
            elif '100MSa' in file:
                source=0.1
                subsample = 1 
                
                realsr=source/subsample
                appendix=f'[{realsr:2.1f} GSa]'
                print("OK - 100 MSa   ",end='')
            else:
                source=0.1
                subsample = 1 
                print(file)
                input('WARNING! NONSTANDARD WF SAMPLING.\nPress key to acknowledge and resume: ')
                realsr=source/subsample
                appendix='[NONSTANDARD SR]'
                print("Warning - Nonstandard SR   ",end='')
      
            
            name = file[:-8]
            print(name,end=' ---> ')
            nameX = file[:-8]+'_{X}.csv'
            nameY = file[:-8]+'_{Y}.csv'
            sourceX = np.loadtxt(nameX)
            print('X ',end='')
            sourceY = np.loadtxt(nameY)
            print('Y ',end='')
            
            crunchedX = sourceX[::subsample].copy()
            crunchedY = sourceY[::subsample].copy()
            ret = np.transpose(np.stack((crunchedX,crunchedY)))
            np.save(name+'__currSR'+appendix+'_{XY}.npy',ret)
            print('DONE')
            if del_csvs:
                os.remove(file[:-8]+'_{X}.csv')
                os.remove(file[:-8]+'_{Y}.csv')

def PlotOverview(set_folder,MA_len=1,axis_upscale=1e9):
    os.chdir(set_folder)
    _, _, filenames = next(os.walk(os.getcwd()), (None, None, []))
    enumr = 0
    enums = 0
    npy_filecount = 0
    
    for file in filenames:
        if file[-4:] == '.npy':
            npy_filecount = npy_filecount+1
    
    rows,cols = factor_int_2(npy_filecount)
    if rows < 4:
        while rows < 4:
            npy_filecount = npy_filecount+1
            rows,cols = factor_int_2(npy_filecount)
    
    fig,ax = plt.subplots(rows,cols,figsize=(16,9))
    for fname in filenames:
        #print(file)
        if fname[-4:] == '.npy':
                        
            if enumr==rows:
                enumr=0
                enums=enums+1
            
            loaded = np.load(fname)
            xx = loaded[:,0]            
            yy = loaded[:,1]
            
            dt = (xx[100]-xx[0])/100
            xxFrom0 = np.arange(0,dt*xx.shape[0],dt)
        
            yyMA = MovingAverage(yy,MA_len)
            xxMA = xx[0:yyMA.shape[0]]             
            xxFrom0 = xxFrom0[0:yyMA.shape[0]] 
        
            ax[enumr,enums].plot(xxFrom0*axis_upscale, yyMA,color='xkcd:crimson',label=fname[0:15])
            ax[enumr,enums].set_xlabel(r't [ns]')
            ax[enumr,enums].set_ylabel(r'OSC readout [V]')
            ax[enumr,enums].legend(loc=1,fontsize='x-small')
            enumr = enumr + 1

    #figname='RL_PulseProp_VariableVCSEL-I_Subs10x'
    #fig.suptitle(figname+f'__MA_{MA_len_t}s__Dev#4-1550nm__16-700-B(-1)',y=0.97)
    fig.tight_layout()
    fig.savefig('FIG__Overview.png')

#%%
if __name__ == '__main__':
    root = Tk()
    root.withdraw()
    root.update()
    set_folder = filedialog.askdirectory(mustexist=True)
    root.destroy()
              
    CrunchCSVs(set_folder)
    #PlotOverview(set_folder,axis_upscale=1e6,MA_len = 1)
#%%
        
#%%   
"""     
fig,ax = plt.subplots(1,1,figsize=(7,4))
workdir = 'D:\\Work\\RTDs - Experimental\\21-09-17_TUe-RTDs_MOD_2'
PlotSingle(os.path.join(workdir,'05o_RL_SAM(1,1)DEVO1C1R3_+0.968V_-6dB_VCS3tunning=2.20mA_(wf1)_[500GSa]__currSR[10.0 GSa]_{XY}.npy'),fig,ax)
"""