# Welcome to StrathLab

This is used to directly control most of our lab instruments, perform measurements in a streamlined way.

The template allows to:

create arbitrary waveforms via custom code,
send waveforms to AWG, controll all key parameters, save those as metadata with each measurement,
readout Keithley state,
readout all oscilloscope channels simultaneously at each shot (via custom R&S readout wrapper)
perform arbitrary number of repetitions over any number of OSC channels,
save and store all the data and metadata in structured, compressed datafiles.
Among other features, the repeated acquisition enables fast and easy recording of large amounts of repeated readouts across multiple oscilloscope channels with all the channels being in sync. Average readout traces are also directly provided for each channel over all acquired repetitions.