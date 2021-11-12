# CLS REIXS Analysis

This is a library to analyse, plot, and export REIXS beamline data. The package is meant to provide a framework to load data into jupyter and enable data interaction.

Further [beamline information](https://reixs.leightsource.ca/) is available on the Website of the Canadian Light Source.

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
sca.load(basedir,'FileName.dat','x_stream','y_stream',1,2,3,4,norm=True)
sca.addScans(basedir,'FileName.dat','x_stream','y_stream',1,2,3,4,norm=False,avg=False)
sca.subtractScans(basedir,'FileName.dat','x_stream','y_stream',1,2,3,4,norm=False,avg=False)
sca.plot()
sca.exporter()
```

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

5. List all scans to analyse (comma-separated)

6. Set optional flags. Options include:

- _norm_ (Normalizes to [0,1])
- _xcoffset_ (Defines a constant shift in the x-stream)
- _xoffset_ (Takes a list of tuples and defines a polynomial fit of the x-stream)
- _ycoffset_ (Defines a constant shift in the y-stream)
- _yoffset_ (Takes a list of tuples and defines a polynomial fit of the y-stream)
  e.g. offset = [(100,102),(110,112),(120,121)]
- _background_ (Subtracts a XEOL background from XEOL scans)
- _energyloss_ (Requires the incident photon energy and moves the resultant MCP scale to energy loss)

#### Absorption Scans

```
xas = XASLoader()
#xas.load(basedir,'Plate2a.dat','TEY',1,4,6,norm=True)
#xas.load(basedir,'Plate2a.dat','PFY[O]',1,4,norm=True)
xas.add(basedir,'Plate2a.dat','PFY[O]',1,4,norm=False,avg=False)
xas.subtract(basedir,'Plate2a.dat','PFY[O]',1,4,6,norm=False,avg=False)
xas.plot()
xas.exporter()
```

#### Emission Scans (MCP)

```
xes = XESLoader()
# Options: XES, rXES
xes.load(basedir,'Plate2a.dat','XES',3,xoffset=[(510,520)])
xes.load(basedir,'Plate2a.dat','XES',3)
xes.load(basedir,'Plate2a.dat','rXES[520:560]',4)
xes.add(basedir,'Plate2a.dat','XES',1,4,)
xes.subtract(basedir,'Plate2a.dat','XES',1,4)
xes.plot()
xes.exporter()
```

#### XRF Scans (SDD)

```
xrf = XRFLoader()
# Options XRF,rXRF
xrf.load(basedir,'Plate2a.dat','XRF',3,xoffset=[(510,520)])
xrf.load(basedir,'Plate2a.dat','XRF',3)
xrf.load(basedir,'Plate2a.dat','rXRF[520:560]',4)
xrf.add(basedir,'Plate2a.dat','XRF',1,4,)
xrf.subtract(basedir,'Plate2a.dat','XRF',1,4)
xrf.plot()
xrf.exporter()
```

#### XEOL Scans

```
xeol = XEOLLoader()
#Options: XEOL, rXEOL
xeol.load(basedir,'RIXS_ES_QA.dat','XEOL',1,2,3,4,background=3)
xeol.plot()
```

### 2d Images

Note: Can only load one scan at a time!

```
load2d = Load2d()
load2d.load(basedir,'Plate2a.dat','Mono Energy','SDD Energy','SDD',1)
load2d.plot()
load2d.exporter()
```

### EEMs (normalized by mesh current)

Note: Can only load one scan at a time!

```
eems = EEMsLoader()
eems.load(basedir,'Plate2a.dat','SDD',1)
eems.load(basedir,'Plate2a.dat','MCP',1)
eems.load(basedir,'RIXS_ES_QA.dat','XEOL',2,background=3)
eems.plot()
eems.exporter()
```

### Mesh Scans

```
mesh = LoadMesh()
mesh.load(basedir,'mesh_scan.txt','Y','Z','TEY',24,norm=True)
mesh.plot()
mesh.exporter()
```

