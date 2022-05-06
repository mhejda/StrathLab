# -*- coding: utf-8 -*-
"""
Library of utility functions for measurements in Antonio Hurtados lab at IoP, Strathclyde
v.01
@author: Matej Hejda, Dafydd Owen-Newns
"""

import warnings
import numpy as np
import os
#import pandas as pd
from time import sleep
from time import time
import matplotlib.pyplot as plt
from IPython.display import Markdown #for text coloring
#import random
#import pickle

try:
    import pyarbtools
except:
    pass

'''
try:
    import pyvisa as visa
except:
    print("Dependency pyvisa not loaded for LabControlToolkit.")       
'''

##############################################################################

##### UTILITY
def Find_Nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def Filesize_Fmt(num, suffix="B"):
    # credit: https://stackoverflow.com/questions/1094841/get-human-readable-version-of-file-size
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"

def Get_Files_From_Dir(directory=None):
    # Returns all filenames in a directory
    if directory == None:
        _, _, filenames = next(os.walk(os.getcwd())) 
        return filenames
    else:
        _, _, filenames = next(os.walk(os.path.join(os.getcwd(),directory)))
        return filenames        

def MovingAverage(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def Clamp(n, minval, maxval):
    try:
        for i, value in enumerate(n):
            n[i]=max(minval, min(maxval, value))
    except:
        n = max(min(n, maxval), minval)
    return n

def Factor_Int_2(n):
    #source: https://stackoverflow.com/questions/39248245/factor-an-integer-to-something-as-close-to-a-square-as-possible
    val = ceil(sqrt(n))
    while True:
        if not n%val:
            val2 = n//val
            break
        val -= 1
    return val, val2

#### SPECIFIC USE UTILITY

def Visualise_Data_Dict(d,lvl=0):
    # by binnev, from: https://stackoverflow.com/questions/15023333/simple-tool-library-to-visualize-huge-python-dict
    # go through the dictionary alphabetically. Modified for StrathLab measurement-data dictionaries
    for k in sorted(d):

        # print the table header if we're at the beginning
        if lvl == 0 and k == sorted(d)[0]:
            print('{:<25} {:<15} {:<10}'.format('KEY','LEVEL','TYPE'))
            print('-'*79)

        indent = '  '*lvl # indent the table to visualise hierarchy
        
        tname = type(d[k]).__name__
        t = f'<{tname}>'

        printable_types = ('int','float','str','tuple')
        if tname in printable_types:
            if 'xpar' in k:
                outp = t+f": {d[k][2]/(d[k][1]-d[k][0]):.2e} Sa/s, {d[k][2]} Sa"
                print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,outp))    
            else:
                print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t+": "+str(d[k])))


        else:
            if k =='readout_osc' or k =='measurement':
                try:
                    try:
                        xpar = d['readout_osc']['mean']['xpar']
                    except:   
                        xpar = d['measurement']['mean']['xpar']
                    rep = d['repeats']
                    outp = t+f": {rep}x | {xpar[2]/(xpar[1]-xpar[0]):.2e} Sa/s, {xpar[2]} Sa"
                    print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,outp))
                except:
                    print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t+": "+str(d[k])))
            else:
                # print details of each entry
                print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t))
            

        # if the entry is a dictionary
        if type(d[k])==dict and k != 'readout_osc' and k != 'measurement':
            # visualise THAT dictionary with +1 indent
            Visualise_Data_Dict(d[k],lvl+1)

#%% MOCK INSTRUMENT (for )
def Initiate_Mock(address='0'):
    mock = Mock()
    return mock

class Mock():
    def __init__(self):
        pass
    def Acq_OSC_Trace(rth, chan,verbose=False):
        return None
    def Keithley_getReadout(sm):
        return None

#%% KEITHLEY 2450 control

def Initiate_Keithley(address="USB0::0x05E6::0x2450::04497230::INSTR"):
    from pymeasure.instruments.keithley import Keithley2450
    
    sm = Keithley2450(address)
    print('Connected to: ',end='')
    print(sm.id)
    return sm

def Keithley_GetReadout(sm):
    ret = sm.ask(':READ? "defbuffer1", READ, SOUR')
    ret = ret.split(',')
    I = eval(ret[0])
    V = eval(ret[1])
    return (I,V)
        
#%% oscilloscope functions

##### ROHDE SCHWARZ OSCILLOSCOPE
### TODO: change name to init/initialise
def Initiate_OSC(address='USB0::0x0AAD::0x0197::1320.5007k08-100963::INSTR'):
    # Creates oscilloscope instance
    from RsInstrument.RsInstrument import RsInstrument
    rth=None
    try:
        # adjust the VISA Resource string to fit your instrument
        rth = RsInstrument(address, True, False)
        #rth = RsInstrument(instr_list[selected], True, False)
        # rth = RsInstrument('USB0::0x0AAD::0x012F::1317.5000K02/103176::INSTR', True, False)
        rth.visa_timeout = 20000 # Timeout for VISA Read Operations
        rth.opc_timeout = 3000 # Timeout for opc-synchronised operations
        rth.instrument_status_checking = True # Error check after each command
        print("Connected to oscilloscope.")
    except Exception as ex:
        print('Error initializing the instrument session:\n' + ex.args[0])
        #exit()
    
    print(f'Device IDN: {rth.idn_string}')
    
    #rth.write_str('ACQuire:SRR 20E9')
    print(f'Device Options: {",".join(rth.instrument_options)}\n')  
    
    rth.data_chunk_size = 100000
    return rth

def Set_OSC_Timebase(rth, left, right):
    dt=right-left
    t0=0.5*(right+left)
    rth.write_str(f"TIM:RANG {dt}")
    rth.write_str(f"TIMebase:HORizontal:POSition {t0}")

def Acq_OSC_Trace(rth, chan,verbose=False,parametric_x=False):
    # Acquires a single trace from a given channel of the scope.
    # Returned format is either linspace-compatible tuple of three values (parametric_x = True),
    #  or the full vector of t-values (parametric_x = False)
    if chan in (1,2,3,4):
        tvalues = rth.query_str(f'CHAN{chan}:WAV1:DATA:HEAD?')
        tvalues = [float(s) for s in tvalues.split(',')]

        scpi_string = f'FORM REAL,32;:CHAN{chan}:DATA?'
        if verbose:
            print(f'Processing SCPI comm: {scpi_string}')
        data_bin_ch = rth.query_bin_or_ascii_float_list(scpi_string)   
        if verbose:
            print(f'OSC: Channel {chan} readout successful.')
        
        if parametric_x:
            return (tvalues[0],tvalues[1],int(tvalues[2])), data_bin_ch
        else:
            t=np.arange(tvalues[0],tvalues[1],(tvalues[1]-tvalues[0])/tvalues[2])
            return t, data_bin_ch
    else:
        raise Warning('Wrong channel selected for readout. Allowed values: 1,2,3,4')
        

#%% awg functions
        
def Initialise_AWG(address="130.159.91.79"):
    awg=None
    try:
        awg = pyarbtools.instruments.M8190A(address, port=5025, timeout=10, reset=False)
        print(awg.instId)
    except:
        Warning("Connecting to AWG failed.")
    return awg

def Calculate_WF_xpar(wf,sr):
    return (0,(wf.shape[0]-1)/sr,wf.shape[0])

def Send_WFs_to_AWG(wf1,
                    out1,
                    ch1_V, 
                    wf2, 
                    out2, 
                    ch2_V, 
                    isHot, 
                    awg, sr, 
                    channels, 
                    rounding_fact=1):
    
    wf1_y = wf1
    wf2_y = wf2
    ch1_voltage = ch1_V
    ch2_voltage = ch2_V
    samplerate = sr
    using_AWG_channels = channels
    wf_rounding_factor = rounding_fact
    
    if using_AWG_channels == (1,2):
        wf1len = wf1_y.shape[0]
        wf2len = wf2_y.shape[0]
        if (wf1len-wf1len%wf_rounding_factor) > (wf2len-wf2len%wf_rounding_factor):
            wf1_y = wf1_y[0:(wf2len-wf2len%wf_rounding_factor)]
            wf2_y = wf2_y[0:(wf2len-wf2len%wf_rounding_factor)]
            print('- Waveform lengths matched: WF1 shortened to match WF2.')
        else:
            wf1_y = wf1_y[0:(wf1len-wf1len%wf_rounding_factor)]
            wf2_y = wf2_y[0:(wf1len-wf1len%wf_rounding_factor)]
            print('- Waveform lengths matched: WF2 shortened to match WF1.')
        print(f'Waveform WF1 length after rounding (samples): {wf1_y.shape[0]}')
        print(f'Waveform WF2 length after rounding (samples): {wf2_y.shape[0]}')
        
        wf1_xpar = Calculate_WF_xpar(wf1_y,samplerate)
        wf2_xpar = Calculate_WF_xpar(wf2_y,samplerate)

    elif using_AWG_channels[0] == 1:
        wf1len = wf1_y.shape[0]
        wf1_y = wf1_y[0:(wf1len-wf1len%wf_rounding_factor)]
        print(f'Waveform WF1 length after rounding (samples): {wf1_y.shape[0]}')
        
        wf1_xpar = Calculate_WF_xpar(wf1_y,samplerate)
        
    elif using_AWG_channels[0] == 2:
        wf2len = wf2_y.shape[0]
        wf2_y = wf2_y[0:(wf2len-wf2len%wf_rounding_factor)]
        print(f'Waveform WF2 length after rounding (samples): {wf2_y.shape[0]}')

        wf2_xpar = Calculate_WF_xpar(wf2_y,samplerate)   
        
    try:
        awg.clear_all_wfm()
    except:
        pass
    
    if isHot:
        awg.configure(fs=samplerate, out1=out1, out2=out2, amp1=ch1_voltage, amp2=ch2_voltage, func1='arb')
        print(f'Sample rate is {awg.fs/1000000} MSa/sec.')
    
        if 1 in using_AWG_channels:
            segment1 = awg.download_wfm(wfmData=wf1_y, ch=1, name='wfm', wfmFormat='real', sampleMkr=1, sampleMkrLength=40)
            awg.play(wfmID=segment1, ch=1)
        if 2 in using_AWG_channels:
            segment2 = awg.download_wfm(wfmData=wf2_y, ch=2, name='wfm', wfmFormat='real', sampleMkr=1, sampleMkrLength=40)
            awg.play(wfmID=segment2, ch=2)
        
        if using_AWG_channels == (1,2):
            fig,ax=plt.subplots(2,1,figsize=(12,5),sharex=True)
    
            ax[0].plot(np.linspace(*wf1_xpar),wf1_y,color='xkcd:indigo',alpha=0.6)
            ax[0].set_facecolor((1.0, 0.47, 0.42,0.2))
            ax[0].set_title('AWG Channel 1 output')    
    
            ax[1].plot(np.linspace(*wf2_xpar),wf2_y,color='xkcd:indigo',alpha=0.6) 
            ax[1].set_title('AWG Channel 2 output')
            ax[1].set_facecolor((1.0, 0.47, 0.42,0.2))
            
            return wf1_xpar, wf2_xpar, fig
        
        elif len(using_AWG_channels) == 1:
            fig,ax=plt.subplots(1,1,figsize=(12,3),sharex=True)
            wf_y = eval(f'wf{using_AWG_channels[0]}_y')
            wf_xpar = eval(f'wf{using_AWG_channels[0]}_xpar')
            
            ax.plot(np.linspace(*wf_xpar),wf_y,color='xkcd:indigo',alpha=0.6) 
            ax.set_title(f'AWG Channel {using_AWG_channels[0]} output')
            ax.set_facecolor((1.0, 0.47, 0.42,0.2))    
            ax.set_xlabel('Time')
            
            if using_AWG_channels[0] == 1:
                return wf1_xpar, None, fig
            else:
                return None, wf2_xpar, fig
    else:
        display (Markdown('<span style="color: #e17701;font-weight:bold">Waveform displayed, but not sent to instrument. Reason: [isHot = False].</span>'))
        if using_AWG_channels == (1,2):
            fig,ax=plt.subplots(2,1,figsize=(12,5),sharex=True)
    
            ax[0].plot(np.linspace(*wf1_xpar),wf1_y,color='xkcd:indigo',alpha=0.6)
            ax[0].set_title('AWG Channel 1 output')    
    
            ax[1].plot(np.linspace(*wf2_xpar),wf2_y,color='xkcd:indigo',alpha=0.6) 
            ax[1].set_title('AWG Channel 2 output')
            
            return wf1_xpar, wf2_xpar, fig
        
        elif len(using_AWG_channels) == 1:
            fig,ax=plt.subplots(1,1,figsize=(12,3),sharex=True)
            
            wf_y = eval(f'wf{using_AWG_channels[0]}_y')
            wf_x = eval(f'wf{using_AWG_channels[0]}_x')
            ax.plot(wf_x,wf_y,color='xkcd:indigo',alpha=0.6) 
            ax.set_title(f'AWG Channel {using_AWG_channels[0]} output')
            ax.set_xlabel('Time')

            if using_AWG_channels[0] == 1:
                return wf1_xpar, None, fig
            else:
                return None, wf2_xpar, fig    
