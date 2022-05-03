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
    warnings.warn("Dependency pyarbtools not loaded for StrathLabToolkit.",RuntimeWarning)
    print()

'''
try:
    import pyvisa as visa
except:
    print("Dependency pyvisa not loaded for LabControlToolkit.")       
'''

##############################################################################

##### UTILITY
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def filesize_fmt(num, suffix="B"):
    # credit: https://stackoverflow.com/questions/1094841/get-human-readable-version-of-file-size
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"

def GetFilesFromDir(directory=None):
    if directory == None:
        _, _, filenames = next(os.walk(os.getcwd())) 
        return filenames
    else:
        _, _, filenames = next(os.walk(os.path.join(os.getcwd(),directory)))
        return filenames        

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

def Keithley_getReadout(sm):
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

def Acq_OSC_Trace(rth, chan,verbose=False):
    if chan in (1,2,3,4):
        tvalues = rth.query_str(f'CHAN{chan}:WAV1:DATA:HEAD?')
        tvalues = [float(s) for s in tvalues.split(',')]
        t=np.arange(tvalues[0],tvalues[1],(tvalues[1]-tvalues[0])/tvalues[2])
        
        scpi_string = f'FORM REAL,32;:CHAN{chan}:DATA?'
        if verbose:
            print(f'Processing SCPI comm: {scpi_string}')
        data_bin_ch = rth.query_bin_or_ascii_float_list(scpi_string)   
        if verbose:
            print(f'OSC: Channel {chan} readout successful.')
        return t, data_bin_ch
    else:
        raise Warning('Wrong channel selected for readout. Allowed values: 1,2,3,4')
        

#%% awg functions
        
def Initialise_AWG(address="130.159.91.79"):
    awg=None
    try:
        awg = pyarbtools.instruments.M8190A(address, port=5025, timeout=10, reset=False)
    except:
        print("could not connect to awg")
    return awg

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
    
    wf1len = wf1_y.shape[0]
    wf2len = wf2_y.shape[0]
    if using_AWG_channels == (1,2):
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
        
        wf1_x = np.arange(0,wf1_y.shape[0]*1/samplerate,1/samplerate)
        wf1_x = wf1_x[0:wf1_y.shape[0]] 
        
        wf2_x = np.arange(0,wf2_y.shape[0]*1/samplerate,1/samplerate)
        wf2_x = wf2_x[0:wf2_y.shape[0]]
    elif using_AWG_channels[0] == 1:
        wf1_y = wf1_y[0:(wf1len-wf1len%wf_rounding_factor)]
        print(f'Waveform WF1 length after rounding (samples): {wf1_y.shape[0]}')
        
        wf1_x = np.arange(0,wf1_y.shape[0]*1/samplerate,1/samplerate)
        wf1_x = wf1_x[0:wf1_y.shape[0]]     
        
    elif using_AWG_channels[0] == 2:
        wf2_y = wf2_y[0:(wf2len-wf2len%wf_rounding_factor)]
        print(f'Waveform WF2 length after rounding (samples): {wf2_y.shape[0]}')
        wf2_x = np.arange(0,wf2_y.shape[0]*1/samplerate,1/samplerate)
        wf2_x = wf2_x[0:wf2_y.shape[0]]    
        
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
    
            ax[0].plot(wf1_x,wf1_y,color='xkcd:indigo',alpha=0.6)
            ax[0].set_facecolor((1.0, 0.47, 0.42,0.2))
            ax[0].set_title('AWG Channel 1 output')    
    
            ax[1].plot(wf2_x,wf2_y,color='xkcd:indigo',alpha=0.6) 
            ax[1].set_title('AWG Channel 2 output')
            ax[1].set_facecolor((1.0, 0.47, 0.42,0.2))
            
            return wf1_x, wf2_x, fig
        
        elif len(using_AWG_channels) == 1:
            fig,ax=plt.subplots(1,1,figsize=(12,3),sharex=True)
            wf_y = eval(f'wf{using_AWG_channels[0]}_y')
            wf_x = eval(f'wf{using_AWG_channels[0]}_x')
            
            ax.plot(wf_x,wf_y,color='xkcd:indigo',alpha=0.6) 
            ax.set_title(f'AWG Channel {using_AWG_channels[0]} output')
            ax.set_facecolor((1.0, 0.47, 0.42,0.2))    
            ax.set_xlabel('Time')
            
            if using_AWG_channels[0] == 1:
                return wf1_x, None, fig
            else:
                return None, wf2_x, fig
    else:
        display (Markdown('<span style="color: #e17701;font-weight:bold">Waveform displayed, but not sent to instrument. Reason: [isHot = False].</span>'))
        if using_AWG_channels == (1,2):
            fig,ax=plt.subplots(2,1,figsize=(12,5),sharex=True)
    
            ax[0].plot(wf1_x,wf1_y,color='xkcd:indigo',alpha=0.6)
            ax[0].set_title('AWG Channel 1 output')    
    
            ax[1].plot(wf2_x,wf2_y,color='xkcd:indigo',alpha=0.6) 
            ax[1].set_title('AWG Channel 2 output')
            
            return wf1_x, wf2_x, fig
        
        elif len(using_AWG_channels) == 1:
            fig,ax=plt.subplots(1,1,figsize=(12,3),sharex=True)
            
            wf_y = eval(f'wf{using_AWG_channels[0]}_y')
            wf_x = eval(f'wf{using_AWG_channels[0]}_x')
            ax.plot(wf_x,wf_y,color='xkcd:indigo',alpha=0.6) 
            ax.set_title(f'AWG Channel {using_AWG_channels[0]} output')
            ax.set_xlabel('Time')

            if using_AWG_channels[0] == 1:
                return wf1_x, None, fig
            else:
                return None, wf2_x, fig    
