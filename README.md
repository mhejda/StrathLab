# StrathLab
Strathclyde Neuromorphic Photonics Lab Control: scripts and notebooks.

# Jupyter Notebook template for data acquisition

Directly control most of lab instruments, create waveforms, send to AWG, readout multiple oscilloscope readouts and save all the data and metadata into structured, compressed datafiles.

## Measurement datafile structure

The measurements are saved in datafiles designed to store all data AND metadata for the measurement.
The goal is to have the datafile completely self-sustained, minimizing the need for referencing the corresponding lab-book entry or writing down parameters and filenames by hand.
The datafile is stored as a LZMA-compressed, serialized (pickled) dictionary object (dict) with following hierarchy:

- `fname` (stores the reference filename defined during acquisition, in case the file would get renamed)
- `date` (full date of acquistion, including HH:MM:SS)
- `repeats` (how many readouts are taken in the measurement)
- `modulation` (stores everything related to AWG modulation)
    - all variables in the scope that include 'ch1_' and/or 'ch_2' are saved alongside the measurement
    - `wf1/2_xpar`
    - `wf1/2_y`
- readout_osc (oscilloscope readout traces, numbered as simple integers)
    - `0`, `1`, `2`, `3`, ...
        - `xpar `
        - `y`
    - `mean` (simple arithmetic mean for all acquired measurements, useful for quickly checking readouts from a datafile)
        - `xpar`
        - `y`

NOTE: `xpar` is a three-value tuple which can be directly fed to `np.linspace(*xpar)` to create a correspodning time-values vector for any recorded data.

Contents of any data-dictionary can be directly viewer with `lab.Visualise_Data_Dict(dict)`.
