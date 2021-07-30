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
sca.exporter('name')
```

1. Specify the variable for the base directory (basedir)
2. Enter the file name of the scan to analyse ('FileName.dat')

3. Options for **x_stream** quantities include:
  * All quantities in the header file
  * *Mono Energy* for the excitation energy
  * *MCP Energy* (uncalibrated)
  * *SDD Energy* (uncalibrated)
  * *XEOL Energy* (uncalibrated, actually the wavelength scale)
  * *Points* (by index)

4. Options for **y_stream** quantities include:
  * All quantities in the header file
  * *TEY* (Total Electron Yield: sample normalized by mesh)
  * *TFY* (Total Fluorescence Yield, normalized by mesh)
  * *PFY* and *iPFY* (Partial Fluorescence Yield and Inverse Partial Fluorescence Yield, both normalized by mesh)
    Specify ROI with brackets, either by XAS edge or energy:
    e.g. *PFY[O]* for PFY at the O K edge
    e.g. *PFY[490:560]* for PFY from 490eV to 560eV
  * *specPFY* (spectrometer PFY, normalized by mesh)
    specify energy range
    e.g. specPFY[500:520]
  * *XES* and *rXES* (X-Ray emission and resonant x-ray emission at selected energies from the spectrometer MCP data)
    e.g. rXES[560:565]
  * *XRF* and *rXRF* (X-Ray fluorescence and resonant x-ray fluorescence at selected energies from the SDD data)
    e.g. rXRF[550:570]
  * *XEOL* and *rXEOL* (XEOL data from the optical spectrometer)
  * *POY* and *TOY* (Partial optical yield and total optical yield, normalized by mesh)
    e.g. POY[300:750]

5. List all scans to analyse (comma-separated)

6. Set optional flags. Options include:
  * *norm* (Normalizes to [0,1])
  * *coffset* (Defines a constant shift in the x-stream)
  * *offset* (Takes a list of tuples and defines a polynomial fit of the x-stream)
    e.g. offset = [(100,102),(110,112),(120,121)]
  * *background* (Subtracts a XEOL background from XEOL scans)

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

#### Emission Scans
```
xes = XESLoader()
```

#### XEOL Scans
```
xeol = XEOLLoader()
```
### 2d Images
To be documented.

### Mesh Scans
To be documented.