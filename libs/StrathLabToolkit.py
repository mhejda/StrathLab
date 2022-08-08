# -*- coding: utf-8 -*-
"""
Library of utility functions for measurements in Antonio Hurtados lab at IoP, Strathclyde
v1.1
@author: Matěj Hejda, Dafydd Owen-Newns
"""
import os, sys, warnings, lzma, pickle

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon

from IPython.display import Markdown #for text coloring
from datetime import datetime
from tqdm import tqdm

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
#### Dictionary subclass, to allow for more convenient access.
#### source: https://goodcode.io/articles/python-dict-object/
class objdict_root(dict):
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)
    
    def plot(self):
        
        fig,ax = plt.subplots(self.readout_osc_count,1,figsize=(6,3+2*self.readout_osc_count))
        #fig.suptitle("Mean traces for all recorded oscilloscope channels.")
        for ind in range(self.readout_osc_count):
            xpar = self[f'readout_osc_{ind}'].xpar
            y = self[f'readout_osc_{ind}'].ymean*1000
            ax[ind].plot(np.linspace(*xpar),y,lw=0.75,color='xkcd:indigo',label=f'readout_osc_{ind}')
            ax[ind].set_ylabel('Voltage [mV]')
            ax[ind].set_xlabel('Time [s]')
            ax[ind].set_xlim(xpar[0],xpar[1])
            ax[ind].legend(fontsize=10)
        
        
    def contents(self,lvl=0,key=None):
        # by binnev, from: https://stackoverflow.com/questions/15023333/simple-tool-library-to-visualize-huge-python-dict
        # go through the dictionary alphabetically. Modified for StrathLab measurement-data dictionaries
        if key==None:
            d=self
        else:
            d=self[key]
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
                if ('readout_osc' in k) or k =='measurement':
                    try:
                        xpar = d[k]['xpar']
                        rep = d['repeats']
                        outp = t+f": {rep}x | {xpar[2]/(xpar[1]-xpar[0]):.2e} Sa/s, {xpar[2]} Sa"
                        print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,outp))
                    except:
                        print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t+": "+str(d[k])))
                else:
                    # print details of each entry
                    print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t))
                
    
            # if the entry is a dictionary
            if (type(d[k])==dict or type(d[k]).__name__=='objdict') and ('readout_osc' not in k) and k != 'measurement':
                # visualise THAT dictionary with +1 indent
                self.contents(lvl=lvl+1,key=k)
            
class objdict(dict):
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)

##### STYLE
"""
import matplotlib.style as style
style.use('default')
plt.rcParams['font.family']='sans-serif'
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['axes.labelsize'] = 11
"""

def Generate_Colormap(colorlist,granularity=1024):
    #credit: @armatita, from https://stackoverflow.com/questions/57268627/matplotlib-color-gradient-between-two-colors
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list(
            "Custom", colorlist, N=granularity)
    return cm

def Adjust_Colormap(cmap, value=1.):
    import colorsys
    import matplotlib.colors as mcolors
    colors = cmap(np.arange(cmap.N))
    hls = np.array([colorsys.rgb_to_hls(*c) for c in colors[:,:3]])
    hls[:,1] *= value
    rgb = np.clip(np.array([colorsys.hls_to_rgb(*c) for c in hls]), 0,1)
    return mcolors.LinearSegmentedColormap.from_list("", rgb)


def Plot_Gradient(ax_object,curvex, curvey,color='xkcd:magenta',ymin=0,ymax=0.02):
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
       
##### UTILITY
def Find_Nearest(array, value):
    # credit: @Demitri at https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
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

def Moving_Average(x, N):
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

def Savefig(fig,fname,path='figs'):
    fig.savefig(f'{path}/FIG__{fname}.pdf')
    fig.savefig(f'{path}/FIG__{fname}.png',dpi=300,transparent=False)

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
            if ('readout_osc' in k) or k =='measurement':
                try:
                    xpar = d[k]['xpar']
                    rep = d['repeats']
                    outp = t+f": {rep}x | {xpar[2]/(xpar[1]-xpar[0]):.2e} Sa/s, {xpar[2]} Sa"
                    print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,outp))
                except:
                    print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t+": "+str(d[k])))
            else:
                # print details of each entry
                print("{:<25} {:<15} {:<10}".format(indent+str(k),lvl,t))
            

        # if the entry is a dictionary
        if (type(d[k])==dict or type(d[k]).__name__=='objdict') and ('readout_osc' not in k) and k != 'measurement':
            # visualise THAT dictionary with +1 indent
            Visualise_Data_Dict(d[k],lvl+1)

#%% MOCK INSTRUMENT (for )
def Initialise_Mock(address='0'):
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

def Initialise_Keithley(address="USB0::0x05E6::0x2450::04497230::INSTR"):
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
def Initialise_OSC(address='USB0::0x0AAD::0x0197::1320.5007k08-100963::INSTR'):
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
        #rth.write_bool(f'EXP:WAV:MULT ',False)
        #rth.write_bool(f'EXP:WAV:INCX ',False)
        
        print("Connected to oscilloscope.")
    except Exception as ex:
        print('Error initializing the instrument session:\n' + ex.args[0])
        #exit()
    
    print(f'Device IDN: {rth.idn_string}')
    
    #rth.write_str('ACQuire:SRR 20E9')
    print(f'Device Options: {",".join(rth.instrument_options)}\n')  
    
    rth.data_chunk_size = 100000
    return rth

def Set_OSC_Channels(rth, acq_channels,verbose=False):
    if not set(acq_channels).issubset((1,2,3,4)):
        warnings.warn("At least one of the specified OSC channel numbers is out of allowed range.")
        sys.exit()  
    
    if len(acq_channels) == 1:
        rth.write_bool(f'EXP:WAV:MULT ',False)
        
        if verbose:
            print('Single channel acq. set on OSC.', end='')    
    elif len(acq_channels) > 1:
        rth.write_bool(f'EXP:WAV:MULT ',True)
        if verbose:
            print('Multichannel acq. set on OSC:  ', end='')
        
        #This is done automatically for all active channels   
        for ii in np.arange(1,5):
            if ii in acq_channels:
                rth.write_bool(f'CHAN{ii}:EXP ',True)
                if verbose:
                    print(f'CH{ii} | ', end='')
            #else:
                #rth.write_bool(f'CHAN{ii}:EXP ',False)

    print()

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
           
        rth.write_str(f'EXP:WAV:SOUR C{chan}W1')       
        tvalues = rth.query_str(f'CHAN{chan}:WAV1:DATA:HEAD?')
        tvalues = [float(s) for s in tvalues.split(',')]
        print(tvalues)

        if verbose:
            print(f'Processing SCPI comm: FORM REAL,32;:CHAN{chan}:WAV1:DATA?')
        data_bin_ch = rth.query_bin_or_ascii_float_list(f'FORM REAL,32;:CHAN{chan}:WAV1:DATA?')   
        if verbose:
            print(f'OSC: Channel {chan} readout successful.')
        
        if parametric_x:
            if tvalues[3] == 1:
                return (tvalues[0],tvalues[1],int(tvalues[2])), data_bin_ch                
            else:
                raise Warning('CH:DATA:HEAD?:4 -> Number of values per sample interval !=1. Increasing xpar length.')
                return (tvalues[0],tvalues[1],int(tvalues[2]*tvalues[3])), data_bin_ch 
        else:
            t=np.arange(tvalues[0],tvalues[1],(tvalues[1]-tvalues[0])/tvalues[2])
            return t, data_bin_ch
    else:
        raise Warning('Wrong channel selected for readout. Allowed values: 1,2,3,4')

def Acq_OSC_Traces(rth, acq_channels, verbose=False):
    # UPDATED VERSION of Acq_OSC_Trace, allows for simultaneous multichannel acquisiton
    # Acquires a single trace from a given channel of the scope.
    # parametric_x = is True by default
    
    traces = {}
    tvalues = rth.query_str(f'CHAN{acq_channels[0]}:DATA:HEAD?')
    tvalues = [float(s) for s in tvalues.split(',')]

    for readch in acq_channels:
        data_bin_ch = rth.query_bin_or_ascii_float_list(f'FORM REAL,32;:CHAN{readch}:DATA?')
        if len(data_bin_ch) >= (tvalues[2]-1):
            break
    
    ch_count = len(acq_channels)
    
    if verbose:
        if (int(tvalues[2])*ch_count != len(data_bin_ch)):
            raise Warning(f'OSC readout check: {int(tvalues[2])}*{ch_count} != {len(data_bin_ch)}')
        else:
            print(f'OSC readout check: {int(tvalues[2])}*{ch_count} == {len(data_bin_ch)}. OK')
    
    for chno, iii in enumerate(acq_channels):
        traces[str(iii)] = data_bin_ch[chno::ch_count]
    return (tvalues[0],tvalues[1],int(tvalues[2])), traces                      

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

def Initialise_Notebook(isHot,datadir='data', figsdir = 'figs'):
    datadir_full = os.path.join(os.getcwd(),datadir)
    # Ensure folders for data and figs exist
    if not os.path.isdir(datadir_full):
        os.mkdir(datadir_full)
    if not os.path.isdir(os.path.join(os.getcwd(),figsdir)):
        os.mkdir(os.path.join(os.getcwd(),figsdir))
    
    # Instrument libraries
    try:
        #from RsInstrument import BinFloatFormat
        #import RsInstrument as rs
        #import pyvisa as visa
        import pyarbtools
        display (Markdown('<span style="color: #1fa774">All instrument control libraries found and loaded.</span>'))
    except:
        display (Markdown('<span style="color: #e17701;font-weight:bold">At least one of the instrument library failed to load. Remote control of lab equipment is not operational.</span>'))

    
    # Initiate remote control for lab instruments
    awg = None
    rth = None
    sm = None
    
    if isHot:
        try:
            awg = pyarbtools.instruments.M8190A('130.159.91.79', port=5025, timeout=30, reset=False)
            print('Connected to AWG.')
            print(awg.query('*IDN?'))
            print()
        except:
            try:
                awg = pyarbtools.instruments.M8190A('127.0.0.1', port=5025, timeout=30, reset=False)
                print('Connected to AWG.')
                print(awg.query('*IDN?'))
                print()
            except:
                warnings.warn('AWG not connected.',RuntimeWarning)
                awg = Mock()
        
        try:
            rth = Initialise_OSC()
        except:
            warnings.warn('Oscilloscope not connected.',RuntimeWarning)
            rth = Mock()
            
        try:
            sm = Initialise_Keithley()
            sm.measure_current(current=0.02, auto_range=True)
            sm.apply_voltage(voltage_range=2, compliance_current=0.02)    
        except:
            warnings.warn('Keithley not connected.',RuntimeWarning)
            sm = Mock()       
            
    else:
        awg = Mock()
        rth = Mock()
        sm = Mock()
    return awg, rth, sm, datadir_full

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
            wf_xpar = eval(f'wf{using_AWG_channels[0]}_xpar')
            ax.plot(np.linspace(*wf_xpar),wf_y,color='xkcd:indigo',alpha=0.6) 
            ax.set_title(f'AWG Channel {using_AWG_channels[0]} output')
            ax.set_xlabel('Time')

            if using_AWG_channels[0] == 1:
                return wf1_xpar, None, fig
            else:
                return None, wf2_xpar, fig    

if __name__ == "__main__":
    rth = Initialise_OSC()
    acq_channels = (3,4)
    Set_OSC_Channels(rth, acq_channels)
    xx,xy = Acq_OSC_Traces(rth, acq_channels,verbose = True)
    plt.plot(np.linspace(*xx),xy[str(acq_channels[0])])

def Get_Modulation_Variables(globals_in,
                             using_AWG_channels,
                             samplerate,
                             saved_obj):
    #### Saving all modulation parameters via keywords search
    if using_AWG_channels != None:
        saved_obj['modulation'] = objdict()
        saved_obj['modulation']['samplerate'] = samplerate
        for jjj in using_AWG_channels:
            #wholex = eval(f'wf{jjj}_x')
            saved_obj['modulation'][f'wf{jjj}_xpar'] = globals_in[f'wf{jjj}_xpar']
            saved_obj['modulation'][f'wf{jjj}_y'] = globals_in[f'wf{jjj}_y']

            print(f'Mod parameters (CH{jjj}) saved:\n[',end='')
            #saved_wf_par_counter = 0
            for var in globals_in:
                if f'ch{jjj}_' in var:
                    saved_obj['modulation'][var] = globals_in[var]
                    #saved_wf_par_counter = saved_wf_par_counter+1
                    print(var,end=', ')
            #print(f'], in total {saved_wf_par_counter} variables saved')
    else:
        print('Modulation is off. No mod parameters saved.')
        saved_obj['modulation'] = None
    return saved_obj
  
def Get_OSC_readouts(acq_channels,
                     repeats,
                     osc_channel_info,
                     saved_obj,
                     awg, 
                     rth, 
                     sm,
                     fname,
                     datadir_full,
                     isSaving):
    #### Filename synthesis
    now = datetime.now()
    saved_obj['date'] = now.strftime("%Y/%m/%d, %H:%M:%S")
    
    full_fname = now.strftime("%Y-%m-%d")+'__'+fname+'.pkl.lz'
    full_filepath = os.path.join(datadir_full,full_fname)
    
    saved_obj['fname'] = fname
    saved_obj['repeats'] = repeats 
    
    #### Catches integer-only channel values, converts them to required tuple format
    if isinstance(acq_channels,int):
        channel_number = acq_channels
        acq_channels = (channel_number,)
    channelcount = len(acq_channels)
    saved_obj['readout_osc_count'] = channelcount
    
    #### Overwrite check
    if os.path.isfile(full_filepath):
        warnings.warn("File exists already! Press Enter to continue...")
        entered = input()
        if entered != '':
            sys.exit()  
        print('Proceeding with acquisiton...')
    
    display (Markdown(f'<span style="color: #014d4e;font-weight:bold">Recording data into: {fname}.</span>'))
    
    #### Input channels check:
    if not set(acq_channels).issubset((1,2,3,4)):
        warnings.warn("At least one of the specified OSC channel numbers is out of allowed range.")
        sys.exit() 
    else:
        Set_OSC_Channels(rth, acq_channels, verbose = True)    
        
    #### Saving Keithley voltage if available
    try:
        sm_I,sm_V = Keithley_GetReadout(sm)
        saved_obj['params']['V_sm'] = sm_V
        print(f'Keithley SM voltage: {sm_V}')
    except:
        pass
    
    #### Create plots: overlay of all measurements, mean+min/max
    if channelcount == 1:
        figR,axR = plt.subplots(2,1,figsize=(9,8),gridspec_kw={'height_ratios':(1,1.5)})
    else:
        figR,axR = plt.subplots(2,channelcount,figsize=(7+2*channelcount,8),gridspec_kw={'height_ratios':(1,1.5)})
 
    #### Get readout shape to create according arrays
    ytotal = objdict()
    ymin = {}
    ymax = {}

    t, ys = Acq_OSC_Traces(rth, acq_channels)
    for ch_no, ch in enumerate(acq_channels):        
        y = np.asarray(ys[str(ch)])

        ytotal[str(ch)] = np.zeros(y.shape)
        ymin[str(ch)] = y.copy()
        ymax[str(ch)] = y.copy() 
        saved_obj[f'readout_osc_{ch_no}'] = objdict()
        try:
            saved_obj[f'readout_osc_{ch_no}']['description'] = osc_channel_info[str(ch)]
        except:
            pass
    
    #### Measure
    pbar = tqdm(total=repeats*len(acq_channels))
    #pbar.set_description()
    for jjj in range(repeats):          
        xpar, ys = Acq_OSC_Traces(rth, acq_channels)
        for ch_no, ch in enumerate(acq_channels):  
            y = np.asarray(y)
            
            saved_obj[f'readout_osc_{ch_no}'][f'xpar'] = xpar
            saved_obj[f'readout_osc_{ch_no}'][f'y{jjj}'] = np.asarray(ys[str(ch)])
            
            if channelcount == 1:
                axR[0].plot(np.linspace(*xpar),ys[str(ch)],color='xkcd:black',lw=2,alpha=1/repeats)
            else:
                axR[0,ch_no].plot(np.linspace(*xpar),ys[str(ch)],color='xkcd:black',lw=2,alpha=1/repeats)
     
            ytotal[str(ch)] = ytotal[str(ch)]+ys[str(ch)]
            ymin[str(ch)] = np.minimum(ymin[str(ch)],ys[str(ch)])
            ymax[str(ch)] = np.maximum(ymax[str(ch)],ys[str(ch)])
            pbar.update(1)
    pbar.close()       

    #### Render mean, maximum, minimum readouts
    
    if len(acq_channels) == 1:
        ytotal[str(ch)] = ytotal[str(ch)]/repeats     
        #axR[1].plot(t,,color='xkcd:evergreen',lw=0.5,alpha=0.5,label='Min')
        axR[1].fill_between(np.linspace(*t),ymin[str(ch)],ymax[str(ch)],color='xkcd:ocean blue',alpha=0.5,label='Min-Max')
        axR[1].plot(np.linspace(*t),ytotal[str(ch)],color='xkcd:black',lw=1,label='Mean trace')
        axR[1].legend()
        axR[0].set_facecolor((1.0, 0.47, 0.42,0.2))
        axR[1].set_facecolor((1.0, 0.47, 0.42,0.2))
        axR[0].set_ylabel('Overlay of all recorded traces')
        axR[1].set_ylabel('Average trace\n + sample-wise min/max')

        #### Dump mean trace into the file for convenience
        #saved_obj[f'readout_osc_{ch_no}']['mean'] = {}
        saved_obj[f'readout_osc_{ch_no}']['ymean'] = ytotal[str(ch)]
    else:    
        for ch_no, ch in enumerate(acq_channels):       
            ytotal[str(ch)] = ytotal[str(ch)]/repeats     
            #axR[1].plot(t,,color='xkcd:evergreen',lw=0.5,alpha=0.5,label='Min')
            axR[1,ch_no].fill_between(np.linspace(*t),ymin[str(ch)],ymax[str(ch)],color='xkcd:ocean blue',alpha=0.5,label='Min-Max')
            axR[1,ch_no].plot(np.linspace(*t),ytotal[str(ch)],color='xkcd:black',lw=1,label='Mean trace')
            axR[1,ch_no].legend()
            axR[0,ch_no].set_facecolor((1.0, 0.47, 0.42,0.2))
            axR[1,ch_no].set_facecolor((1.0, 0.47, 0.42,0.2))
            axR[0,ch_no].set_ylabel('Overlay of all recorded traces')
            axR[1,ch_no].set_ylabel('Average trace\n + sample-wise min/max')

            #### Dump mean trace into the file for convenience
            #saved_obj[f'readout_osc_{ch_no}']['mean'] = {}
            #saved_obj[f'readout_osc_{ch_no}']['mean']['xpar'] = xpar
            saved_obj[f'readout_osc_{ch_no}']['ymean'] = ytotal[str(ch)]
        figR.tight_layout()

    #### Save everything into LZMA-compressed pickled (serialized) file
    if isSaving:
        #memsize = sys.getsizeof(saved_obj)
        #print(f'Object size (memory): {memsize} [{lab.filesize_fmt(memsize)}]')
        with lzma.open(full_filepath,"wb",preset=3) as f:
            pickle.dump(saved_obj,f)
            print(f'Measurement saved, time: {saved_obj["date"]}')
        print('Filesize: '+Filesize_Fmt(os.stat(full_filepath).st_size)) 
        
        figR.savefig(full_filepath[:-7]+".png")
    
    return saved_obj