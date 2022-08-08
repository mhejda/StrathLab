# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:12:13 2021

@author: thisi
"""

import numpy as np
import matplotlib.pyplot as plt
import os, sys, importlib, lzma, pickle

sys.path.insert(0,'C:/git/StrathLab/libs') #Lab PC
sys.path.insert(0,'D:/Repositories/StrathLab/libs') #XMG
sys.path.insert(0,'D:/Repositories/PyRTD/libs') #XMG

import StrathLabToolkit as lab
importlib.reload(lab)

def Convert(remove_csvs = True):
    files = lab.Get_Files_From_Dir()
    
    converted = 0
    
    for fname in files:
        if '{Y}' in fname and os.path.exists(fname[:-7]+'{Xpar}.csv'):
            with lzma.open((os.path.join(os.getcwd(),fname[:-14]+'.pkl.lz')),"wb",preset=3) as f:
                d = lab.objdict()
                
                print('Found (X,Y): '+ fname[:-7],end='')
                
                jjj = 0
                d['fname'] = fname
                d[f'repeats'] = 1
                d[f'modulation'] = None
                d[f'readout_osc_0'] = {}
                xpar_arr = np.loadtxt(fname[:-7]+'{Xpar}.csv',delimiter=",")
                d[f'readout_osc_0'][f'xpar'] = (xpar_arr[0],xpar_arr[1],int(xpar_arr[2]))
                d[f'readout_osc_0'][f'y{jjj}'] = np.loadtxt(fname,delimiter=",")
                
                pickle.dump(d,f)
                
                #plt.plot(np.linspace(*d[f'readout_osc_0'][f'xpar']),d[f'readout_osc_0'][f'y{jjj}'])
            print(' --- OK')
            
            if remove_csvs:
                os.remove(fname)
                os.remove(fname[:-7]+'{Xpar}.csv')
            converted = converted + 1
    print(f'= Converted {converted} X+Y pairs. =')
#%%

if __name__ == '__main__':
    Convert()
    