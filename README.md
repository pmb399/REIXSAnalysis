# CLS REIXS Analysis

This is a library to analyse, plot, and export REIXS beamline data. The package is meant to provide a framework to load data into jupyter and enable data interaction.

Further [beamline information](https://reixs.leightsource.ca/) is available on the Website of the Canadian Light Source.

## Installation
Install the package from PyPi with the pip package manager. This is the recommended way to obtain a copy for your local machine and will install all required dependencies.
```
    $ pip install reixs
```
You will also need [Jupyter Notebook](https://github.com/jupyter) together with python 3 on your local machine.

In the case that certain widgets aren't rendered properly, make sure yo enable the appropriate jupyter extensions
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

## Setup necesarry inputs
from reixs.LoadData import *
from bokeh.io import show, output_notebook
output_notebook(hide_banner=True)
```

### 1d plots
#### General Loader1d
```
sca = Load1d()
sca.load(...)
sca.addScans(...)
sca.subtractScans(...)
sca.plot()
sca.exporter(...)
```

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