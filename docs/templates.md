# Operation & Templates

This is intended to be used with Jupyter Notebook, two templates are provided for this.

* _Acquisition template contains code for sending data to the AWG and reading from the oscilloscope
* _Processing template can be used to view data previously collected

## Acquisition

The cells in this notebook are used for sending waveforms to the AWG and then reading from the oscilloscope. On reading, all information used (waveforms, device settings and user defined data) is saved alongside the collected data.

### AWG waveforms

Here the parameters of the AWG can be set, and data can be sent as waveforms.

The channels to activate are specified in a tuple, if only using a single channel an int can be used.

Regions where user data enter the procedure are marked with comments.


### Acquisition

This defines the data to be recorded and other data to be saved.

The object saved_obj is a dictionary that will be serialised and saved, this will automatically contain all information relevant to the recording it is produced for. Any extra information can be added as needed.

