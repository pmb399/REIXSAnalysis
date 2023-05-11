# CLS REIXS Analysis

This is a library to analyse, plot, and export REIXS beamline data. The package is meant to provide a framework to load data into jupyter and enable data interaction.

Further [beamline information](https://reixs.lightsource.ca/) is available on the Website of the Canadian Light Source.

## Installation

Install the package from PyPi with the pip package manager. This is the recommended way to obtain a copy for your local machine and will install all required dependencies.

```
    $ pip install reixs
```

You will also need [Jupyter Notebook](https://github.com/jupyter) together with python 3 on your local machine.

In case that certain widgets aren't rendered properly, make sure to enable the appropriate jupyter extensions

```
    $ jupyter nbextension enable --py widgetsnbextension
```

## Running

Launch your local jupyter installation with

```
    $ jupyter notebook
```

## Examples

### Load the required module

Before you start, you will need to import the required reixs package, enable bokeh plotting, and set the base directory.

```
## Define base directory
basedir = "/home/braun/ownCloud/Beamtime/example_data/"

## Setup necessarry inputs
from reixs.LoadData import *
from bokeh.io import show, output_notebook
output_notebook(hide_banner=True)
```

### 1d plots

#### General Loader1d

```
sca = Load1d()
sca.load(basedir,'FileName.dat','x_stream','y_stream',1,2,3,4)  # Loads multiple scans individually
sca.add(basedir,'FileName.dat','x_stream','y_stream',1,2,3,4)  # Adds multiple scans
sca.subtract(basedir,'FileName.dat','x_stream','y_stream',1,2,3,4,norm=False) # Subtracts scans from the first scan
sca.xlim(lower_lim,upper_lim) # Sets the horizontal axis plot region
sca.ylim(lower_lim,upper_lim) # Sets the vertical axis plot region
sca.plot_legend("pos string as per bokeh") # Determines a specific legend position
sca.vline(position) # Draws a vertical line
sca.hline(position) # Draws a horizontal line
sca.label(pos_x,pos_y,'Text') # Adds a label to the plot
sca.plot() # Plots the defined object
sca.exporter() # Exports the data by calling an exporter widget
```

0. Create "Loader" object

1. Specify the variable for the base directory (basedir)

2. Enter the file name of the scan to analyse ('FileName.dat')

3. Options for **x_stream** quantities include:
- All quantities in the header file
- _Mono Energy_ for the excitation energy
- _MCP Energy_ (uncalibrated)
- _SDD Energy_ (uncalibrated)
- _XEOL Energy_ (uncalibrated, actually the wavelength scale)
- _Points_ (by index)

4. Options for **y_stream** quantities include:
- All quantities in the header file
- _TEY_ (Total Electron Yield: sample normalized by mesh)
- _TFY_ (Total Fluorescence Yield, normalized by mesh)
- _PFY_ and _iPFY_ (Partial Fluorescence Yield and Inverse Partial Fluorescence Yield, both normalized by mesh)
  Specify ROI with brackets, either by XAS edge or energy:
  e.g. _PFY[O]_ for PFY at the O K edge
  e.g. _PFY[490:560]_ for PFY from 490eV to 560eV
- _specPFY_ (spectrometer PFY, normalized by mesh)
  specify energy range
  e.g. specPFY[500:520]
- _XES_ and _rXES_ (X-Ray emission and resonant x-ray emission at selected energies from the spectrometer MCP data)
  e.g. rXES[560:565]
- _XRF_ and _rXRF_ (X-Ray fluorescence and resonant x-ray fluorescence at selected energies from the SDD data)
  e.g. rXRF[550:570]
- _XEOL_ and _rXEOL_ (XEOL data from the optical spectrometer)
- _POY_ and _TOY_ (Partial optical yield and total optical yield, normalized by mesh)
  e.g. POY[300:750]
- _EY_ or _Sample_ (Sample current, not normalized by mesh)
- _Mesh_ (Mesh current)
- _ET_ (Energy Transfer data, integrates over energy loss ROI and probes constant final states, sometimes referred to as CET scan)
  specify energy transfer region
  e.g. ET[-2:5] to probe mostly scattering close to the elastic line
- _rLOSS_ (Resonantly excited emission data on energy loss scale, integrates over incident energy ROIS and probes constant intermediate states, sometimes referred to as CIE scan)
  specify incident energy region
  e.g. rLOSS[620:640]

5. List all scans to analyse (comma-separated)

6. Set optional flags. Options include:
- _norm_ (Normalizes to [0,1])
- _xcoffset_ (Defines a constant shift in the x-stream)
- _xoffset_ (Takes a list of tuples and defines a polynomial fit of the x-stream)
- _ycoffset_ (Defines a constant shift in the y-stream)
- _yoffset_ (Takes a list of tuples and defines a polynomial fit of the y-stream)
  e.g. offset = [(100,102),(110,112),(120,121)]
- _background_ (Subtracts a XEOL background from XEOL scans)
  Set to True, uses the getXEOLback function with the background data stored (only supported with HDF5)
  Specify scan number, subtracts the XEOL scan taken at this particular scan
- _energyloss_ (Transfers the resultant MCP scale to energy loss 
  Set to True, then takes mean of mono energy array
  Specify float with the incident photon energy
- _grid_x_ (Takes a list with three arguments to apply 1d interpolation gridding)
  e.g. grid_x = [Start Energy, Stop Energy, Delta]
- _savgol_ (Takes a list with two or three arguments to apply data smoothing and derivatives)
  e.g. savgol = [Window length, Polynomial order, deriavtive] as specified in the scipy Savitzky-Golay filter
- _binsize_ (int, allows to perform data binning to improve Signal-to-Noise)
- _legend_items_ (dict={scan_number:"name"}, overwrites generic legend names; works for the _load_ method)
- _legend_item_ (str, overwrites generic legend name in the _add_/_subtract_ method)

#### Absorption Scans

```
xas = XASLoader()
xas.load(basedir,'Plate2a.dat','TEY',1,4,6)
xas.load(basedir,'Plate2a.dat','PFY[O]',1,4)
xas.add(basedir,'Plate2a.dat','PFY[500:520]',1,4)
xas.subtract(basedir,'Plate2a.dat','specPFY[500:520]',1,4,6)
xas.plot()
xas.exporter()
```

#### Emission Scans (MCP)

```
xes = XESLoader()
# Options: XES, rXES
xes.load(basedir,'Plate2a.dat','XES',3,xoffset=[(515,520),(520,525),(530,535)])
xes.load(basedir,'Plate2a.dat','rXES[520:560]',4)
xes.add(basedir,'Plate2a.dat','XES',1,4)
xes.subtract(basedir,'Plate2a.dat','XES',1,4)
xes.plot()
xes.exporter()
```

#### XRF Scans (SDD)

```
xrf = XRFLoader()
# Options XRF,rXRF
xrf.load(basedir,'Plate2a.dat','XRF',3)
xrf.load(basedir,'Plate2a.dat','rXRF[520:560]',4)
xrf.add(basedir,'Plate2a.dat','XRF',1,4,)
xrf.subtract(basedir,'Plate2a.dat','XRF',1,4)
xrf.plot()
xrf.exporter()
```

#### XEOL Scans (Optical Spectrometer)

```
xeol = XEOLLoader()
#Options: XEOL, rXEOL
xeol.load(basedir,'RIXS_ES_QA.dat','XEOL',1,2,3,4,background=3)
xeol.load(basedir,'RIXS_ES_QA.dat','XEOL',1,2,3,4,background=True)
xeol.plot()
```

### 2d Images

#### General loader for MCA detector data

Note: Can only load one scan at a time!

```
load2d = Load2d()
load2d.load(basedir,'Filename.dat','x_stream','y_stream','detector',1)
load2d.plot()
load2d.exporter()
```

0. Create "Loader" object

1. Specify the variable for the base directory (basedir)

2. Enter the file name of the scan to analyse ('FileName.dat')

3. Options for **x_stream** quantities include:
- All quantities in the header file
- _Mono Energy_ for the excitation energy

4. Options for **y_stream** quantities include:
- _SDD Energy_ (Energy scale of the SDD detector)
- _MCP Energy_ (Energy scale of the MCP detector)
- _XEOL Energy_ (Wavelength scale of the XEOL optical spectrometer)

5. Options for **detector** quantities include:
- _SDD_ (SDD detector MCA)
- _MCP_ (MCP detector MCA)
- _XEOL_ (XEOL optical spectrometer MCA)

6. List all scans to analyse (comma-separated)

7. Set optional flags. Options include:
- _norm_ (Normalizes to [0,1])
- _xcoffset_ (Defines a constant shift in the x-stream)
- _xoffset_ (Takes a list of tuples and defines a polynomial fit of the x-stream)
- _ycoffset_ (Defines a constant shift in the y-stream)
- _yoffset_ (Takes a list of tuples and defines a polynomial fit of the y-stream)
  e.g. offset = [(100,102),(110,112),(120,121)]
- _background_ (Subtracts a XEOL background from XEOL scans)
  Set to True, uses the getXEOLback function with the background data stored (only supported with HDF5)
  Specify scan number, subtracts the XEOL scan taken at this particular scan
- _energyloss_ (Transfers the excitation-emission map to energy loss scale
- _grid_x_ (Takes a list with three arguments to apply 1d interpolation gridding)
  e.g. grid_x = [Start Energy, Stop Energy, Delta]

#### EEMs (normalized by mesh current, special version of the general 2d image loader)

Note: Can only load one scan at a time!

```
eems = EEMsLoader()
eems.load(basedir,'Plate2a.dat','SDD',1)
eems.load(basedir,'Plate2a.dat','MCP',1)
eems.load(basedir,'RIXS_ES_QA.dat','XEOL',2,background=3)
eems.plot()
eems.exporter()
```

### Mesh Scans (Plots a 2d histogram)

```
mesh = LoadMesh()
mesh.load(basedir,'Filename.txt','x_stream','y_stream','z_stream',24)
mesh.plot()
mesh.exporter()
```

0. Create "Loader" object

1. Specify the variable for the base directory (basedir)

2. Enter the file name of the scan to analyse ('FileName.dat')

3. Options for **x_stream** quantities include:
- All quantities in the header file
- _Mono Energy_ for the excitation energy
- _SDD Energy_ (Energy scale of the SDD detector)
- _MCP Energy_ (Energy scale of the MCP detector)
- _XEOL Energy_ (Wavelength scale of the XEOL optical spectrometer)

4. Options for **y_stream** quantities include:
- All quantities in the header file
- _Mono Energy_ for the excitation energy
- _SDD Energy_ (Energy scale of the SDD detector)
- _MCP Energy_ (Energy scale of the MCP detector)
- _XEOL Energy_ (Wavelength scale of the XEOL optical spectrometer)

5. Options for **z_stream** quantities include:
- All quantities in the header file
- All special quantities as specified for the Load1d() function

6. List all scans to analyse (comma-separated)

7. Set optional flags. Options include:
- _norm_ (Normalizes to [0,1])
- _xcoffset_ (Defines a constant shift in the x-stream)
- _xoffset_ (Takes a list of tuples and defines a polynomial fit of the x-stream)
- _ycoffset_ (Defines a constant shift in the y-stream)
- _yoffset_ (Takes a list of tuples and defines a polynomial fit of the y-stream)
  e.g. offset = [(100,102),(110,112),(120,121)]
- _background_ (Subtracts a XEOL background from XEOL scans)
  Set to True, uses the getXEOLback function with the background data stored (only supported with HDF5)
  Specify scan number, subtracts the XEOL scan taken at this particular scan
