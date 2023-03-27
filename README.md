# StrathLab
Strathclyde Neuromorphic Photonics Lab Control: scripts and notebooks.

Full documentation at https://strathlab.readthedocs.io/en/latest/

# Jupyter Notebook template for data acquisition

**Main idea:** Directly control most of our lab instruments, perform measurements in a streamlined way. 

The template allows to:
- create arbitrary waveforms via custom code, 
- send waveforms to AWG, controll all key parameters, save those as metadata with each measurement, 
- readout Keithley state, 
- readout all oscilloscope channels simultaneously at each shot (via custom R&S readout wrapper)
- perform arbitrary number of repetitions over any number of OSC channels,
- save and store all the data and metadata in structured, compressed datafiles. 

Among other features, the repeated acquisition enables fast and easy recording of large amounts of repeated readouts across multiple oscilloscope channels with all the channels being in sync. Average readout traces are also directly provided for each channel over all acquired repetitions.

## Measurement datafile structure

The measurements are saved in datafiles designed to store all data AND metadata for the measurement.
The goal is to have the datafile completely self-sustained, minimizing the need for referencing the corresponding lab-book entry or writing down parameters and filenames by hand.
The datafile is stored as a LZMA-compressed, serialized (pickled), modified dictionary subclass (objdict), which wraps dict with ability to use Class-like queries for stored data (f.e. objdict.repeats). It has the following hierarchy:

- `fname` _(stores the reference filename defined during acquisition, in case the file would get renamed)_
- `date` _(full date of acquistion, including HH:MM:SS)_
- `repeats` _(how many readouts are taken in the measurement)_
- `modulation` _(stores everything related to AWG modulation)_
    - _all variables in the scope that include `ch1_` and/or `ch2_` are saved alongside the measurement_
    - `wf1/2_xpar`
    - `wf1/2_y`
- `readout_osc_{n}` _(oscilloscope readout traces, numbered as simple integers, where {n} is an integer enumerating the acquired channel. For single channel readout, just set to 0)_
    - `xpar `
    - `y0`, `y1`, `y2`, `y3`, ...
    - `ymean` _(simple arithmetic mean for all acquired measurements, useful for quickly checking readouts from a datafile)_

NOTE: `xpar` is a three-value tuple which can be directly fed to `np.linspace(*xpar)` to create a correspodning time-values vector for any recorded data.

Contents of any data-dictionary can be directly visualized and viewed with `lab.Visualise_Data_Dict(dict)`.
