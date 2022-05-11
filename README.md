# StrathLab
Strathclyde Neuromorphic Photonics Lab Control: scripts and notebooks.

# Jupyter Notebook template for data acquisition

Directly control most of lab instruments, create waveforms, send to AWG, readout multiple oscilloscope readouts and save all the data and metadata into structured, compressed datafiles.

## Measurement datafile structure

The measurements are saved in datafiles designed to store all data AND metadata for the measurement.
The goal is to have the datafile completely self-sustained, minimizing the need for referencing the corresponding lab-book entry or writing down parameters and filenames by hand.
The datafile is stored as a LZMA-compressed, serialized (pickled) dictionary object (dict) with following hierarchy:

- `fname` _(stores the reference filename defined during acquisition, in case the file would get renamed)_
- `date` _(full date of acquistion, including HH:MM:SS)_
- `repeats` _(how many readouts are taken in the measurement)_
- `modulation` _(stores everything related to AWG modulation)_
    - _all variables in the scope that include `ch1_` and/or `ch2_` are saved alongside the measurement_
    - `wf1/2_xpar`
    - `wf1/2_y`
- `readout_osc` _(oscilloscope readout traces, numbered as simple integers)_
    - `0`, `1`, `2`, `3`, ...
        - `xpar `
        - `y`
    - `mean` _(simple arithmetic mean for all acquired measurements, useful for quickly checking readouts from a datafile)_
        - `xpar`
        - `y`

NOTE: `xpar` is a three-value tuple which can be directly fed to `np.linspace(*xpar)` to create a correspodning time-values vector for any recorded data.

Contents of any data-dictionary can be directly viewer with `lab.Visualise_Data_Dict(dict)`.