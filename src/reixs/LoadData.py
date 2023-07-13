# Scientific Modules
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, interp2d

# Plotting
from bokeh.io import push_notebook
from bokeh.plotting import show, figure
from bokeh.models import ColumnDataSource, HoverTool, LinearColorMapper, LogColorMapper, ColorBar, Span, Label

# Utilities
import os
from collections import defaultdict
import io
import shutil

# Widgets
import ipywidgets as widgets
from IPython.display import display
from ipyfilechooser import FileChooser

# Edge Dict
from .edges import EdgeDict

# Data Processing Functions
from .util import COLORP, all_list_entries_equal
from .xeol import *
from .add_subtract import ScanAddition, ScanSubtraction, ImageAddition, ImageSubtraction
from .sca import loadSCAscans
from .mca import loadMCAscans
from .mesh import loadMeshScans
from .rsxs_mcp import loadRSXS1dROIscans, loadRSXS2dROIscans, loadRSXSImageStack
from .beamline_info import loadSCAbeamline, get_single_beamline_value

#########################################################################################
#########################################################################################


class Load1d:
    """Class to load generic 1d (x,y) data."""

    def __init__(self):
        self.data = list()
        self.type = list()
        self.x_stream = list()
        self.filename = list()
        self.plot_lim_x = [":", ":"]
        self.plot_lim_y = [":", ":"]
        self.legend_loc = 'outside'
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()

    def load(self, basedir, file, x_stream, y_stream, *args, **kwargs):
        """
        Load one or multiple specific scan(s) for selected streams.

        Parameters
        ----------
        basedir : string
            Specifiy the absolute or relative path to experimental data.
        file : string
            Specify the file name (either ASCII or HDF5).
        x_stream : string
            Specifiy the data for the horizontal axis.
            Use: "Mono Energy", "MCP Energy", "SDD Energy", "XEOL Energy", "Points", or any SCA scalar array.
        y_stream : string
            Specifiy the data for the vertical axis.
            Use: "TEY", "TFY, "PFY", "iPFY", "XES", "rXES", "specPFY",
                  "XRF", "rXRF", "XEOL", "rXEOL", "POY", "TOY", "EY", "Sample", "Mesh", "ET", or any SCA scalar array.
        *args : int
            Separate scan numbers with comma.
        **kwargs: multiple, optional
            Options:
                norm : boolean
                    Norm the spectra to [0,1].
                    default: True
                xoffset : list of tuples
                    Offset the x-axis by applying a polynomial fit.
                    default: None
                xcoffset : float
                    Offset x-axis by constant value.
                    default : None 
                yoffset : list of tuples
                    Offset the y-axis by applying a polynomial fit.
                    default : None 
                ycoffset : float
                    Offset y-axis by constant value.
                    default : None
                background : int or boolean
                    Apply background selection for XEOL. 
                    Select True when using with getXEOLback function or select background scan number.
                    default : None
                energyloss : boolean or float
                    Convert emission energy to energy loss.
                    Select True to extract incident mono energy from file or specify float manually.
                    default : None
                grid_x: list
                    Grid data evenly on specified grid [low,high,step size]
                    default: [None, None, None]
                savgol : tupel
                    Apply a Savitzky-Golay filter for smoothing, and optionally take derivatives.
                    arg1 : savgol window length
                    arg2 : savgol polynomial order
                    arg 3: derivative order, optional
                    default : None, 
                binsize : int
                    Bin data by reducing data points via averaging.
                    Int must be exponent of 2.
                    default : None
        """

        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", True)
        kwargs.setdefault("is_XAS", False)
        # Append all REIXS scan objects to scan list in current object.
        self.data.append(loadSCAscans(
            basedir, file, x_stream, y_stream, *args, **kwargs))
        self.type.append(y_stream)
        self.x_stream.append(x_stream)
        self.filename.append(file)

    def add(self, basedir, file, x_stream, y_stream, *args, **kwargs):
        """
        Add specified scans for selected streams.

        Parameters
        ----------
        See loader function.
        Adds all scans specified in *args.
        """

        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", False)
        kwargs.setdefault("avg", False)
        kwargs.setdefault("is_XAS", False)
        # Append all REIXS scan objects to scan list in current object.
        self.data.append(ScanAddition(
            basedir, file, x_stream, y_stream, *args, **kwargs))
        self.x_stream.append(x_stream)
        self.type.append(y_stream)
        self.filename.append(file)

    def subtract(self, basedir, file, x_stream, y_stream, *args, **kwargs):
        """
        Subtract specified scans for selected streams.

        Parameters
        ----------
        See loader function.
        Subtracts all scans from the first element. May add scans in first element by specifying list of scans as first *arg.

        """
        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", False)
        kwargs.setdefault("is_XAS", False)
        # Append all REIXS scan objects to scan list in current object.
        self.data.append(ScanSubtraction(
            basedir, file, x_stream, y_stream, *args, **kwargs))
        self.x_stream.append(x_stream)
        self.type.append(y_stream)
        self.filename.append(file)

    def xlim(self, lower, upper):
        """
        Set x-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_x[0] = lower
        self.plot_lim_x[1] = upper

    def ylim(self, lower, upper):
        """
        Set y-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_y[0] = lower
        self.plot_lim_y[1] = upper

    def plot_legend(self, pos):
        """
        Overwrite default legend position.

        Parameters
        ----------
        pos : string
            See bokeh manual for available options.
        """
        self.legend_loc = pos

    def vline(self, pos, **kwargs):
        """
        Draw a vertical line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_vlines.append([pos, kwargs])

    def hline(self, pos, **kwargs):
        """
        Draw a horizontal line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_hlines.append([pos, kwargs])

    def label(self, pos_x, pos_y, text, **kwargs):
        """
        Draw a text box in the plot.

        Parameters
        ----------
        pos_x : float
        pos_y : float
        text : string
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_labels.append([pos_x, pos_y, text, kwargs])

    def plot(self, linewidth=4, title=None, xlabel=None, ylabel=None, plot_height=450, plot_width=700, y_axis_type='linear'):
        """
        Plot all data assosciated with class instance/object.

        Parameters
        ----------
        linewidth : int, optional
        title : string, optional
        xlabel : string, optional
        ylabel : string, optional
        plot_height : int, optional
        plot_width : int, optional
        """

        # Organize all data assosciated with object in sorted dictionary.
        plot_data = defaultdict(list)
        for i, val in enumerate(self.data):
            for k, v in val.items():
                plot_data["x_stream"].append(v.x_stream)
                plot_data["y_stream"].append(v.y_stream)
                plot_data['x_name'].append(self.x_stream[i])
                plot_data['y_name'].append(self.type[i])
                plot_data['filename'].append(self.filename[i])
                plot_data['scan'].append(v.scan)
                plot_data['legend'].append(v.legend)

        # Get the colours for the glyphs.
        numlines = len(plot_data['scan'])
        plot_data['color'] = COLORP[0:numlines]

        source = ColumnDataSource(plot_data)

        # Set up the bokeh plot
        p = figure(height=plot_height, width=plot_width,y_axis_type=y_axis_type,
                   tools="pan,wheel_zoom,box_zoom,reset,crosshair,save")
        p.multi_line(xs='x_stream', ys='y_stream', legend_group="legend",
                     line_width=linewidth, line_color='color', line_alpha=0.6,
                     hover_line_color='color', hover_line_alpha=1.0,
                     source=source)

        # Set up the information for hover box
        p.add_tools(HoverTool(show_arrow=False, line_policy='next', tooltips=[
            ('Scan', '@scan'),
            ('File', '@filename'),
            ("(x,y)", "(@x_name, @y_name)"),
            ("(x,y)", "($x, $y)")
        ]))

        p.toolbar.logo = None

        # Overwrite plot properties if requested.
        if self.legend_loc == 'outside':
            p.add_layout(p.legend[0], 'right')
        else:
            p.legend.location = self.legend_loc

        if self.plot_lim_y[0] != ':':
            p.y_range.start = self.plot_lim_y[0]
        if self.plot_lim_y[1] != ':':
            p.y_range.end = self.plot_lim_y[1]

        if self.plot_lim_x[0] != ':':
            p.x_range.start = self.plot_lim_x[0]
        if self.plot_lim_x[1] != ':':
            p.x_range.end = self.plot_lim_x[1]

        if len(self.plot_hlines) > 0:
            for line_props in self.plot_hlines:
                line = Span(location=line_props[0],
                            dimension='width', **line_props[1])
                p.add_layout(line)

        if len(self.plot_vlines) > 0:
            for line_props in self.plot_vlines:
                line = Span(
                    location=line_props[0], dimension='height', **line_props[1])
                p.add_layout(line)

        if len(self.plot_labels) > 0:
            for label_props in self.plot_labels:
                label = Label(
                    x=label_props[0], y=label_props[1], text=label_props[2], **label_props[3])
                p.add_layout(label)

        if title != None:
            p.title.text = str(title)
        if xlabel != None:
            p.xaxis.axis_label = str(xlabel)
        if ylabel != None:
            p.yaxis.axis_label = str(ylabel)
        show(p)

    def get_data(self):
        """Make data available in memory as exported to file.

        Returns
        -------
        dfT : pandas DataFrame 
            All loaded data.
        files : list
            List of all loaded files.
        """
        
        files = list()
        series_data = list()
        series_header = list()

        # Iterate over all "load" calls
        for i, val in enumerate(self.data):
            # Iterate over all scans per load call.
            for k, v in val.items():
                name = f"~{self.filename[i]}"
                if name not in files:
                    files.append(name)
                fileindex = files.index(name)

                # Append the x_stream data and header name
                series_data.append(pd.Series(v.x_stream))
                series_header.append(f"F{fileindex+1}_S{v.scan}_I{i+1}-{self.x_stream[i]}")

                # Append the y_stream data and header name
                series_data.append(pd.Series(v.y_stream))
                series_header.append(f"F{fileindex+1}_S{v.scan}_I{i+1}-{self.type[i]}")

        dfT = pd.DataFrame(series_data).transpose(copy=True)
        dfT.columns = series_header

        return dfT, files

    def export(self, filename):
        """
        Export and write data to specified file.

        Parameters
        ----------
        filename : string
        """
        dfT, files = self.get_data()

        # Open file.
        with open(f"{filename}.csv", 'w') as f:
            string = '# '
            # Generate list of files for legend.
            for idx, file in enumerate(files):
                string += f"F{idx+1} {file},"
            string += '\n'
            f.write(string)
            # Write pandas dataframe to file.
            dfT.to_csv(f, index=False, line_terminator='\n')

        print(f"Successfully wrote DataFrame to {filename}.csv")

    def exporter(self):
        """Interactive exporter widget."""
        current_dir = os.path.dirname(os.path.realpath("__file__"))

        self.exportfile = FileChooser(current_dir)
        self.exportfile.use_dir_icons = True
        #self.exportfile.filter_pattern = '*.csv'

        button = widgets.Button(
            description='Save data file',
            disabled=False,
            button_style='info',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='Save data to file',
            icon='save'  # (FontAwesome names without the `fa-` prefix)
        )

        button.on_click(self.exportWidgetStep)
        display(self.exportfile, button)

    def exportWidgetStep(self, my):
        # Helper function for exporter widget.
        file = os.path.join(self.exportfile.selected_path,
                            self.exportfile.selected_filename)
        self.export(file)

# The following classes all inherit from Load1d and provide quick-access to frequently used spectroscopy.


class XASLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = 'Mono Energy'
        kwargs.setdefault("norm", True)
        kwargs["is_XAS"] = True
        super().load(basedir, file, x_stream, y_stream, *args, **kwargs)

    def plot(self, **kwargs):
        kwargs.setdefault(
            "title", 'Absorption spectra normalized by mesh current')
        kwargs.setdefault("xlabel", "Incident Photon Energy (eV)")
        super().plot(**kwargs)

    def add(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "Mono Energy"
        kwargs.setdefault("avg", True)
        kwargs.setdefault("norm", True)
        kwargs["is_XAS"] = True
        super().add(basedir, file, x_stream, y_stream, *args, **kwargs)

    def subtract(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "Mono Energy"
        kwargs.setdefault("norm", True)
        kwargs["is_XAS"] = True
        super().subtract(basedir, file, x_stream, y_stream, *args, **kwargs)


class XESLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, **kwargs):
        kwargs.setdefault("norm", False)
        x_stream = "MCP Energy"
        super().load(basedir, file, x_stream, y_stream, *
                     args, **kwargs)

    def plot(self, **kwargs):
        kwargs.setdefault("title", 'Summed MCP emission spectra')
        kwargs.setdefault("xlabel", "Uncalibrated MCP Energy (eV)")
        kwargs.setdefault("ylabel", "Counts (arb. units)")
        super().plot(**kwargs)

    def add(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "MCP Energy"
        kwargs.setdefault("avg", False)
        kwargs.setdefault("norm", False)
        super().add(basedir, file, x_stream, y_stream, *args, **kwargs)

    def subtract(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "MCP Energy"
        kwargs.setdefault("norm", False)
        super().subtract(basedir, file, x_stream, y_stream, *
                         args, **kwargs)


class XRFLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "SDD Energy"
        kwargs.setdefault("norm", False)
        super().load(basedir, file, x_stream, y_stream, *
                     args, **kwargs)

    def plot(self, **kwargs):
        kwargs.setdefault("title", 'Summed SDD emission spectra')
        kwargs.setdefault("xlabel", "Uncalibrated SDD Energy (eV)")
        kwargs.setdefault("ylabel", "Counts (arb. units)")
        super().plot(**kwargs)

    def add(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "SDD Energy"
        kwargs.setdefault("norm", False)
        kwargs.setdefault("avg", False)
        super().add(basedir, file, x_stream, y_stream, *args, **kwargs)

    def subtract(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "SDD Energy"
        kwargs.setdefault("norm", False)
        super().subtract(basedir, file, x_stream, y_stream, *
                         args, **kwargs)


class XEOLLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "XEOL Energy"
        kwargs.setdefault("norm", False)
        super().load(basedir, file, x_stream, y_stream, *args, **kwargs)

    def plot(self, **kwargs):
        kwargs.setdefault("title", 'Summed XEOL spectra')
        kwargs.setdefault("xlabel", "Uncalibrated XEOL Wavelength (nm)")
        kwargs.setdefault("ylabel", "Counts (arb. units)")
        super().plot(**kwargs)

    def add(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "XEOL Energy"
        kwargs.setdefault("norm", False)
        kwargs.setdefault("avg", False)
        super().add(basedir, file, x_stream, y_stream, *args, **kwargs)

    def subtract(self, basedir, file, y_stream, *args, **kwargs):
        x_stream = "XEOL Energy"
        kwargs.setdefault("norm", False)
        super().subtract(basedir, file, x_stream, y_stream, *args, **kwargs)


#########################################################################################

class Load2d:
    """Class to load generic 2d (x,y, z) image data of a detector."""

    def __init__(self):
        self.data = list()
        self.x_stream = list()
        self.y_stream = list()
        self.detector = list()
        self.filename = list()
        self.plot_lim_x = [":", ":"]
        self.plot_lim_y = [":", ":"]
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()

    def load(self, basedir, file, x_stream, y_stream, detector, *args, **kwargs):
        """
        Load one or multiple specific scan(s) for selected streams.

        Parameters
        ----------
        basedir : string
            Specifiy the absolute or relative path to experimental data.
        file : string
            Specify the file name (either ASCII or HDF5).
        x_stream : string
            Specifiy the data for the horizontal axis.
            Use: "Mono Energy" or any SCA scalar array.
        y_stream : string
            Specifiy the data for the vertical axis.
            Use: "MCP Energy", "SDD Energy", "XEOL Energy"
        detector : string
            Use: "MCP", "SDD", or "XEOL".
        *args : int
            Specify one (1) scan to load.
        **kwargs: multiple, optional
            Options:
                norm : boolean
                    Norm the spectra to [0,1].
                    default: True
                xoffset : list of tuples
                    Offset the x-axis by applying a polynomial fit.
                    default: None
                xcoffset : float
                    Offset x-axis by constant value.
                    default : None 
                yoffset : list of tuples
                    Offset the y-axis by applying a polynomial fit.
                    default : None 
                ycoffset : float
                    Offset y-axis by constant value.
                    default : None
                background : int or boolean
                    Apply background selection for XEOL. 
                    Select True when using with getXEOLback function or select background scan number.
                    default : None
                energyloss : boolean or float
                    Convert emission energy to energy loss.
                    Select True to extract incident mono energy from file or specify float manually.
                    default : None
                grid_x: list
                    Grid data evenly on specified grid [low,high,step size]
                    default: [None, None, None]
                grid_y: list
                    Grid data evenly on specified grid [low,high,step size]
                    default: [None, None, None]
        """

        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", False)
        kwargs.setdefault("xoffset", None)
        kwargs.setdefault("xcoffset", None)
        kwargs.setdefault("yoffset", None)
        kwargs.setdefault("ycoffset", None)
        kwargs.setdefault("background", None)
        kwargs.setdefault("grid_x", [None, None, None])
        kwargs.setdefault("grid_y", [None, None, None])
        kwargs.setdefault("energyloss", False)

        # Ensure that only one scan is loaded.
        if len(args) != 1:
            raise TypeError("You may only select one scan at a time")
        if self.data != []:
            raise TypeError("You can only append one scan per object")
        self.data.append(loadMCAscans(basedir, file, x_stream,
                         y_stream, detector, *args, **kwargs))
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.detector.append(detector)
        self.filename.append(file)

        if kwargs['energyloss'] == True:
            self.x_stream[-1] = "Energy loss (eV)"
            self.y_stream[-1] = "Mono Energy (eV)"

    def add(self, basedir, file, x_stream, y_stream, detector, *args, **kwargs):
        """
        Add specified images for selected streams.

        Parameters
        ----------
        See loader function.
        Adds all scans specified in *args.
        """

        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", False)
        kwargs.setdefault("xoffset", None)
        kwargs.setdefault("xcoffset", None)
        kwargs.setdefault("yoffset", None)
        kwargs.setdefault("ycoffset", None)
        kwargs.setdefault("background", None)
        kwargs.setdefault("grid_x", [None, None, None])
        kwargs.setdefault("grid_y", [None, None, None])
        kwargs.setdefault("energyloss", False)

        self.data.append(ImageAddition(basedir, file, x_stream,
                         y_stream, detector, *args, **kwargs))
        
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.detector.append(detector)
        self.filename.append(file)

        if kwargs['energyloss'] == True:
            self.x_stream[-1] = "Energy loss (eV)"
            self.y_stream[-1] = "Mono Energy (eV)"

    def subtract(self, basedir, file, x_stream, y_stream, detector, *args, **kwargs):
        """
        Subtract specified images for selected streams.

        Parameters
        ----------
        See loader function.
        Subtracts all imnages from the first element.

        """
        # Set the defaults if not specified in **kwargs.
        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", False)
        kwargs.setdefault("xoffset", None)
        kwargs.setdefault("xcoffset", None)
        kwargs.setdefault("yoffset", None)
        kwargs.setdefault("ycoffset", None)
        kwargs.setdefault("background", None)
        kwargs.setdefault("grid_x", [None, None, None])
        kwargs.setdefault("grid_y", [None, None, None])
        kwargs.setdefault("energyloss", False)

        # Append all REIXS scan objects to scan list in current object.
        self.data.append(ImageSubtraction(basedir, file, x_stream,
                         y_stream, detector, *args, **kwargs))
        
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.detector.append(detector)
        self.filename.append(file)

        if kwargs['energyloss'] == True:
            self.x_stream[-1] = "Energy loss (eV)"
            self.y_stream[-1] = "Mono Energy (eV)"

    def xlim(self, lower, upper):
        """
        Set x-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_x[0] = lower
        self.plot_lim_x[1] = upper

    def ylim(self, lower, upper):
        """
        Set y-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_y[0] = lower
        self.plot_lim_y[1] = upper

    def vline(self, pos, **kwargs):
        """
        Draw a vertical line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_vlines.append([pos, kwargs])

    def hline(self, pos, **kwargs):
        """
        Draw a horizontal line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_hlines.append([pos, kwargs])

    def label(self, pos_x, pos_y, text, **kwargs):
        """
        Draw a text box in the plot.

        Parameters
        ----------
        pos_x : float
        pos_y : float
        text : string
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_labels.append([pos_x, pos_y, text, kwargs])

    def plot(self, title=None, kind='Image', xlabel=None, ylabel=None, plot_height=600, plot_width=600, 
            vmin=None, vmax=None, colormap = "linear"):
        """
        Plot all data assosciated with class instance/object.

        Parameters
        ----------
        title : string, optional
        kind : string, optional
        xlabel : string, optional
        ylabel : string, optional
        plot_height : int, optional
        plot_width : int, optional
        vmin : float, optional
        vmax : float, optional
        colormap : string
            Use: "linear" or "log"
        """
        # Iterate over the one (1) scan in object - this is for legacy reason and shall be removed in the future.
        for i, val in enumerate(self.data):
            for k, v in val.items():

                # Create the figure
                p = figure(height=plot_height, width=plot_width, tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                           tools="pan,wheel_zoom,box_zoom,reset,hover,crosshair,save")
                p.x_range.range_padding = p.y_range.range_padding = 0

                # Gridded scales now calculated directly during the MCA load and only need to be referenced here

                # must give a vector of image data for image parameter
                if vmin == None:
                    mapper_low = v.new_z.min()
                else:
                    mapper_low = vmin

                if vmax == None:
                    mapper_high = v.new_z.max()
                else:
                    mapper_high = vmax

                if colormap == "linear":
                    myMapper = LinearColorMapper
                elif colormap == "log":
                    myMapper = LogColorMapper
                else:
                    raise UserWarning("Only 'linear' and 'log' implemented.")

                color_mapper = myMapper(palette="Viridis256",
                                                 low=mapper_low,
                                                 high=mapper_high)

                # Plot image and use limits as given by even grid.
                p.image(image=[v.new_z], x=v.xmin, y=v.ymin, dw=v.xmax-v.xmin,
                        dh=v.ymax-v.ymin, color_mapper=color_mapper, level="image")
                p.grid.grid_line_width = 0.5

                # Defining properties of color mapper
                color_bar = ColorBar(color_mapper=color_mapper,
                                     label_standoff=12,
                                     location=(0, 0),
                                     title='Counts')
                p.add_layout(color_bar, 'right')

                # Overwrite plot properties if selected.
                if self.plot_lim_y[0] != ':':
                    p.y_range.start = self.plot_lim_y[0]
                if self.plot_lim_y[1] != ':':
                    p.y_range.end = self.plot_lim_y[1]

                if self.plot_lim_x[0] != ':':
                    p.x_range.start = self.plot_lim_x[0]
                if self.plot_lim_x[1] != ':':
                    p.x_range.end = self.plot_lim_x[1]

                if len(self.plot_hlines) > 0:
                    for line_props in self.plot_hlines:
                        line = Span(
                            location=line_props[0], dimension='width', **line_props[1])
                        p.add_layout(line)

                if len(self.plot_vlines) > 0:
                    for line_props in self.plot_vlines:
                        line = Span(
                            location=line_props[0], dimension='height', **line_props[1])
                        p.add_layout(line)

                if len(self.plot_labels) > 0:
                    for label_props in self.plot_labels:
                        label = Label(
                            x=label_props[0], y=label_props[1], text=label_props[2], **label_props[3])
                        p.add_layout(label)

            if title != None:
                p.title.text = str(title)
            else:
                p.title.text = f'{self.detector[i]} {kind} for Scan {k}'
            if xlabel != None:
                p.xaxis.axis_label = str(xlabel)
            else:
                p.xaxis.axis_label = str(self.x_stream[i])
            if ylabel != None:
                p.yaxis.axis_label = str(ylabel)
            else:
                p.yaxis.axis_label = f"{self.y_stream[i]}"

            p.toolbar.logo = None

            show(p)

    def get_data(self):
        """Make data available in memory as exported to file.

        Returns
        -------
        f : string.IO object
            Motor and Detector Scales. Pandas Data Series.
            1) Rewind memory with f.seek(0)
            2) Load with pandas.read_csv(f,skiprows=3)
        g : string.IO object
            Actual gridded detector image.
            1) Rewind memory with g.seek(0)
            2) Load with numpy.genfromtxt(g,skip_header=4)
        """
        # Set up the data frame and the two string objects for export
        f = io.StringIO()
        g = io.StringIO()
        series_data = list()
        series_header = list()

        for i, val in enumerate(self.data):
            for k, v in val.items():
                # Gridded scales now calculated directly during the MCA load and only need to be referenced here

                # Start writing string f
                f.write("========================\n")
                f.write(
                    f"F~{self.filename[i]}_S{v.scan}_{self.detector[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                f.write("========================\n")

                # Start writing string g
                g.write("========================\n")
                g.write(
                    f"F~{self.filename[i]}_S{v.scan}_{self.detector[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                g.write("========================\n")

                # Append data to string now.
                # Append x-stream
                series_data.append( pd.Series(v.new_x))
                series_header.append("Motor Scale Gridded")

                # Append y-stream
                series_data.append(pd.Series(v.new_y))
                series_header.append("Detector Scale Gridded")

                dfT = pd.DataFrame(series_data).transpose(copy=True)
                dfT.columns = series_header
                dfT.to_csv(f, index=False, line_terminator='\n')

                g.write("=== Image ===\n")
                np.savetxt(g, v.new_z, fmt="%.9g")

            return f, g

    def export(self, filename):
        """
        Export and write data to specified file.

        Parameters
        ----------
        filename : string
        """
        f, g, = self.get_data()

        # Dump both strings in file.
        # Need to rewind memory location of String.IO to move to beginning.
        # Copy string content to file with shutil.
        with open(f"{filename}.txt_scale", "a") as scales:
            f.seek(0)
            shutil.copyfileobj(f, scales)

        with open(f"{filename}.txt_matrix", "a") as matrix:
            g.seek(0)
            shutil.copyfileobj(g, matrix)

        print(f"Successfully wrote Image data to {filename}.txt")

    def exporter(self):
        """Interactive exporter widget."""
        current_dir = os.path.dirname(os.path.realpath("__file__"))

        self.exportfile = FileChooser(current_dir)
        self.exportfile.use_dir_icons = True
        #self.exportfile.filter_pattern = '*.txt'

        button = widgets.Button(
            description='Save data file',
            disabled=False,
            button_style='info',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='Save data to file',
            icon='save'  # (FontAwesome names without the `fa-` prefix)
        )

        button.on_click(self.exportWidgetStep)
        display(self.exportfile, button)

    def exportWidgetStep(self, my):
        # Helper function for exporter widget
        file = os.path.join(self.exportfile.selected_path,
                            self.exportfile.selected_filename)
        self.export(file)


class EEMsLoader(Load2d):
    """Specific 2d loader for excitation-emission-maps."""

    def get_scale(self, detector):
        if detector == "MCP":
            return "MCP Energy"
        elif detector == "SDD":
            return "SDD Energy"
        elif detector == "XEOL":
            return "XEOL Energy"
        else:
            raise TypeError("Detector not defined.")

    def load(self, basedir, file, detector, *args, **kwargs):
        x_stream = 'Mono Energy'
        y_stream = self.get_scale(detector)

        super().load(basedir, file, x_stream, y_stream, detector, *args, **kwargs)

    def add(self, basedir, file, detector, *args, **kwargs):
        x_stream = "Mono Energy"
        y_stream = self.get_scale(detector)

        super().add(basedir, file, x_stream, y_stream, detector, *args, **kwargs)

    def subtract(self, basedir, file, detector, *args, **kwargs):
        x_stream = "Mono Energy"
        y_stream = self.get_scale(detector)

        super().subtract(basedir, file, x_stream, y_stream, detector, *args, **kwargs)

#########################################################################################


class LoadMesh(Load2d):
    """Class to display (x,y,z) scatter data."""

    def __init__(self):
        self.data = list()
        self.x_stream = list()
        self.y_stream = list()
        self.z_stream = list()
        self.filename = list()
        self.plot_lim_x = [":", ":"]
        self.plot_lim_y = [":", ":"]
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()

        # Use this so we can inherit from Load2d for plotting
        self.detector = self.z_stream

    def load(self, basedir, file, x_stream, y_stream, z_stream, *args, **kwargs):

        # Set the defaults if not specified in **kwargs.
        kwargs.setdefault("norm", False)
        kwargs.setdefault("xoffset", None)
        kwargs.setdefault("xcoffset", None)
        kwargs.setdefault("yoffset", None)
        kwargs.setdefault("ycoffset", None)
        kwargs.setdefault("background", None)
        #kwargs.setdefault("grid_x",[None, None, None])
        #kwargs.setdefault("grid_y",[None, None, None])

        self.data.append(loadMeshScans(basedir, file, x_stream,
                         y_stream, z_stream, *args, **kwargs))
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.z_stream.append(z_stream)
        self.filename.append(file)

    def add(self):
        raise UserWarning("Functionality not yet implemented.")
    
    def subtract(self):
        raise UserWarning("Functionality not yet implemented.")
    
    def plot(self, *args, **kwargs):
        kwargs.setdefault('kind', "Histogram")

        super().plot(*args, **kwargs)

    def get_data(self):
        f = io.StringIO()
        g = io.StringIO()
        for i, val in enumerate(self.data):
            for k, v in val.items():
                # Have the gridded data ready now from loader
                f.write("========================\n")
                f.write(
                    f"F~{self.filename[i]}_S{v.scan}_{self.z_stream[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                f.write("========================\n")

                g.write("========================\n")
                g.write(
                    f"F~{self.filename[i]}_S{v.scan}_{self.z_stream[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                g.write("========================\n")

                f.write("=== x-axis bin edges ===\n")
                np.savetxt(f, v.xedge)
                f.write("=== y-axis bin edges ===\n")
                np.savetxt(f, v.yedge)
                g.write("=== Histogram ===\n")
                np.savetxt(g, v.new_z, fmt="%.9g")
        return f, g

    def export(self, filename):

        f, g, = self.get_data()

        with open(f"{filename}.txt_scale", "a") as scales:
            f.seek(0)
            shutil.copyfileobj(f, scales)

        with open(f"{filename}.txt_matrix", "a") as matrix:
            g.seek(0)
            shutil.copyfileobj(g, matrix)

        print(f"Successfully wrote Histogram data to {filename}.txt")

class MeshLoader(LoadMesh):
    pass

#########################################################################################


class ImageROILoader():
    def __init__(self):
        self.data = list()
        self.filename = list()
        self.norm = list()
        self.plot_lim_x = [":", ":"]
        self.plot_lim_y = [":", ":"]
        self.legend_loc = 'outside'
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()

    def load(self, basedir, file, images, axis, *args, x=[None, None], y=[None, None], norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None):
        self.data.append(loadRSXS1dROIscans(basedir, file, images, axis, *args, x=x, y=y,
                         norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv))
        self.norm.append(norm)
        self.filename.append(file)

    def xlim(self, lower, upper):
        """
        Set x-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_x[0] = lower
        self.plot_lim_x[1] = upper

    def ylim(self, lower, upper):
        """
        Set y-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_y[0] = lower
        self.plot_lim_y[1] = upper

    def plot_legend(self, pos):
        """
        Overwrite default legend position.

        Parameters
        ----------
        pos : string
            See bokeh manual for available options.
        """
        self.legend_loc = pos

    def vline(self, pos, **kwargs):
        """
        Draw a vertical line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_vlines.append([pos, kwargs])

    def hline(self, pos, **kwargs):
        """
        Draw a horizontal line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_hlines.append([pos, kwargs])

    def label(self, pos_x, pos_y, text, **kwargs):
        """
        Draw a text box in the plot.

        Parameters
        ----------
        pos_x : float
        pos_y : float
        text : string
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_labels.append([pos_x, pos_y, text, kwargs])

    def plot(self, linewidth=4, title='Image ROI', xlabel=None, ylabel='Counts (arb. units)', plot_height=450, plot_width=700):
        plot_data = defaultdict(list)
        for i, val in enumerate(self.data):
            for k, v in val.items():
                for image, imagevalue in v.caxis.items():
                    plot_data["x_stream"].append(imagevalue)
                    plot_data["y_stream"].append(v.imagesca[image])
                    plot_data['x_name'].append(v.caxis_label)
                    plot_data['y_name'].append("Counts (arb. units)")
                    plot_data['filename'].append(self.filename[i])
                    plot_data['scan'].append(v.scan)
                    plot_data['image'].append(image)
                    plot_data['legend'].append(f"S{v.scan}_I{image}")

        numlines = len(plot_data['scan'])
        plot_data['color'] = COLORP[0:numlines]

        source = ColumnDataSource(plot_data)

        p = figure(height=plot_height, width=plot_width,
                   tools="pan,wheel_zoom,box_zoom,reset,crosshair,save")
        p.multi_line(xs='x_stream', ys='y_stream', legend_group="legend",
                     line_width=linewidth, line_color='color', line_alpha=0.6,
                     hover_line_color='color', hover_line_alpha=1.0,
                     source=source)

        p.add_tools(HoverTool(show_arrow=False, line_policy='next', tooltips=[
            ('Scan', '@scan'),
            ('Image', '@image'),
            ('File', '@filename'),
            ("(x,y)", "(@x_name, @y_name)"),
            ("(x,y)", "($x, $y)")
        ]))

        # Overwrite plot properties if requested.
        if self.legend_loc == 'outside':
            p.add_layout(p.legend[0], 'right')
        else:
            p.legend.location = self.legend_loc

        if self.plot_lim_y[0] != ':':
            p.y_range.start = self.plot_lim_y[0]
        if self.plot_lim_y[1] != ':':
            p.y_range.end = self.plot_lim_y[1]

        if self.plot_lim_x[0] != ':':
            p.x_range.start = self.plot_lim_x[0]
        if self.plot_lim_x[1] != ':':
            p.x_range.end = self.plot_lim_x[1]

        if len(self.plot_hlines) > 0:
            for line_props in self.plot_hlines:
                line = Span(location=line_props[0],
                            dimension='width', **line_props[1])
                p.add_layout(line)

        if len(self.plot_vlines) > 0:
            for line_props in self.plot_vlines:
                line = Span(
                    location=line_props[0], dimension='height', **line_props[1])
                p.add_layout(line)

        if len(self.plot_labels) > 0:
            for label_props in self.plot_labels:
                label = Label(
                    x=label_props[0], y=label_props[1], text=label_props[2], **label_props[3])
                p.add_layout(label)

        p.toolbar.logo = None

        if all_list_entries_equal(plot_data['x_name']) == True:
            xlabel = plot_data['x_name'][0]
        else:
            xlabel = "Warning: Mixed quantity"

        if title != None:
            p.title.text = str(title)
        if xlabel != None:
            p.xaxis.axis_label = str(xlabel)
        if ylabel != None:
            p.yaxis.axis_label = str(ylabel)
        show(p)

    def get_data(self):
        
        files = list()
        series_data = list()
        series_header = list()

        for i, val in enumerate(self.data):
            for k, v in val.items():
                name = f"~{self.filename[i]}"
                if name not in files:
                    files.append(name)
                fileindex = files.index(name)

                for image, imagevalue in v.caxis.items():
                    # Append the x stream
                    series_data.append(pd.Series(imagevalue))
                    series_header.append(f"F{fileindex+1}_S{v.scan}_I{i+1}_Img{image}-{v.caxis_label}")

                    # Append the y stream
                    series_data.append(pd.Series(v.imagesca[image]))
                    series_header.append(f"F{fileindex+1}_S{v.scan}_I{i+1}_Img{image}_Counts")

        dfT = pd.DataFrame(series_data).transpose(copy=True)
        dfT.columns = series_header

        return dfT, files

    def export(self, filename):

        dfT, files = self.get_data()

        with open(f"{filename}.csv", 'w') as f:
            string = '# '
            for idx, file in enumerate(files):
                string += f"F{idx+1} {file},"
            string += '\n'
            f.write(string)
            dfT.to_csv(f, index=False, line_terminator='\n')

        # dfT.to_csv(f"{filename}.csv",index=False)
        print(f"Successfully wrote DataFrame to {filename}.csv")

    def exporter(self):
        current_dir = os.path.dirname(os.path.realpath("__file__"))

        self.exportfile = FileChooser(current_dir)
        self.exportfile.use_dir_icons = True
        #self.exportfile.filter_pattern = '*.csv'

        button = widgets.Button(
            description='Save data file',
            disabled=False,
            button_style='info',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='Save data to file',
            icon='save'  # (FontAwesome names without the `fa-` prefix)
        )

        button.on_click(self.exportWidgetStep)
        display(self.exportfile, button)

    def exportWidgetStep(self, my):
        file = os.path.join(self.exportfile.selected_path,
                            self.exportfile.selected_filename)
        self.export(file)


class StackROILoader():
    def __init__(self):
        self.data = list()
        self.filename = list()
        self.norm = list()
        self.plot_lim_x = [":", ":"]
        self.plot_lim_y = [":", ":"]
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()

    def load(self, basedir, file, axis, *args, x=[None, None], y=[None, None], norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None):
        self.data.append(loadRSXS2dROIscans(basedir, file, axis, *args, x=x, y=y, norm=norm,
                         xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset))
        self.norm.append(norm)
        self.filename.append(file)

    def xlim(self, lower, upper):
        """
        Set x-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_x[0] = lower
        self.plot_lim_x[1] = upper

    def ylim(self, lower, upper):
        """
        Set y-axis plot window limits.

        Parameters
        ----------
        lower : float
        upper : float
        """
        self.plot_lim_y[0] = lower
        self.plot_lim_y[1] = upper

    def vline(self, pos, **kwargs):
        """
        Draw a vertical line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_vlines.append([pos, kwargs])

    def hline(self, pos, **kwargs):
        """
        Draw a horizontal line in the plot.

        Parameters
        ----------
        pos : float
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_hlines.append([pos, kwargs])

    def label(self, pos_x, pos_y, text, **kwargs):
        """
        Draw a text box in the plot.

        Parameters
        ----------
        pos_x : float
        pos_y : float
        text : string
        **kwargs : dict, optional
            See bokeh manual for available options.
        """
        self.plot_labels.append([pos_x, pos_y, text, kwargs])

    def plot(self, title=None, xlabel=None, ylabel=None, plot_height=600, plot_width=600):
        for i, val in enumerate(self.data):
            for k, v in val.items():

                p = figure(height=plot_height, width=plot_width, tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                           tools="pan,wheel_zoom,box_zoom,reset,hover,crosshair,save")
                p.x_range.range_padding = p.y_range.range_padding = 0

                # must give a vector of image data for image parameter
                color_mapper = LinearColorMapper(palette="Viridis256")

                p.image(image=[v.imagemca], x=v.mcp_x.min(), y=v.mcp_y.min(), dw=v.mcp_x.max()-v.mcp_x.min(),
                        dh=v.mcp_y.max()-v.mcp_y.min(), color_mapper=color_mapper, level="image")
                p.grid.grid_line_width = 0.5

                # Defining properties of color mapper
                color_bar = ColorBar(color_mapper=color_mapper,
                                     label_standoff=12,
                                     location=(0, 0),
                                     title='Counts')
                p.add_layout(color_bar, 'right')

                # Overwrite plot properties if selected.
                if self.plot_lim_y[0] != ':':
                    p.y_range.start = self.plot_lim_y[0]
                if self.plot_lim_y[1] != ':':
                    p.y_range.end = self.plot_lim_y[1]

                if self.plot_lim_x[0] != ':':
                    p.x_range.start = self.plot_lim_x[0]
                if self.plot_lim_x[1] != ':':
                    p.x_range.end = self.plot_lim_x[1]

                if len(self.plot_hlines) > 0:
                    for line_props in self.plot_hlines:
                        line = Span(
                            location=line_props[0], dimension='width', **line_props[1])
                        p.add_layout(line)

                if len(self.plot_vlines) > 0:
                    for line_props in self.plot_vlines:
                        line = Span(
                            location=line_props[0], dimension='height', **line_props[1])
                        p.add_layout(line)

                if len(self.plot_labels) > 0:
                    for label_props in self.plot_labels:
                        label = Label(
                            x=label_props[0], y=label_props[1], text=label_props[2], **label_props[3])
                        p.add_layout(label)

            if title != None:
                p.title.text = str(title)
            else:
                p.title.text = f'Image Stack ROI for Scan {k}'
            if xlabel != None:
                p.xaxis.axis_label = str(xlabel)
            else:
                p.xaxis.axis_label = v.caxis_labels[0]
            if ylabel != None:
                p.yaxis.axis_label = str(ylabel)
            else:
                p.yaxis.axis_label = v.caxis_labels[1]

            p.toolbar.logo = None

            show(p)

    def get_data(self):
        f = io.StringIO()
        g = io.StringIO()

        for i, val in enumerate(self.data):
            for k, v in val.items():
                f.write("========================\n")
                f.write(
                    f"F~{self.filename[i]}_S{v.scan}_{v.caxis_labels[0]}_{v.caxis_labels[1]}\n")
                f.write("========================\n")

                g.write("========================\n")
                g.write(
                    f"F~{self.filename[i]}_S{v.scan}_{v.caxis_labels[0]}_{v.caxis_labels[1]}\n")
                g.write("========================\n")

                f.write("=== Scale x Gridded ===\n")
                np.savetxt(f, v.mcp_x, fmt="%.9g")
                f.write("=== Scale y Gridded ===\n")
                np.savetxt(f, v.mcp_y, fmt="%.9g")
                g.write("=== Image ===\n")
                np.savetxt(g, v.imagemca, fmt="%.9g")

        return f, g

    def export(self, filename):

        f, g, = self.get_data()

        with open(f"{filename}.txt_scale", "a") as scales:
            f.seek(0)
            shutil.copyfileobj(f, scales)

        with open(f"{filename}.txt_matrix", "a") as matrix:
            g.seek(0)
            shutil.copyfileobj(g, matrix)

        print(f"Successfully wrote Image data to {filename}.txt")

    def exporter(self):
        current_dir = os.path.dirname(os.path.realpath("__file__"))

        self.exportfile = FileChooser(current_dir)
        self.exportfile.use_dir_icons = True
        #self.exportfile.filter_pattern = '*.txt'

        button = widgets.Button(
            description='Save data file',
            disabled=False,
            button_style='info',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='Save data to file',
            icon='save'  # (FontAwesome names without the `fa-` prefix)
        )

        button.on_click(self.exportWidgetStep)
        display(self.exportfile, button)

    def exportWidgetStep(self, my):
        file = os.path.join(self.exportfile.selected_path,
                            self.exportfile.selected_filename)
        self.export(file)


class ImageStackLoader():
    def __init__(self):
        self.data = list()
        self.filename = list()
        self.norm = list()

    def load(self, basedir, file, arg):
        if self.data != []:
            raise UserWarning("Can only load one movie at a time.")
        else:
            self.data.append(loadRSXSImageStack(basedir, file, arg))
            self.filename.append(file)

    def plot(self, title=None, xlabel=None, ylabel=None, plot_height=600, plot_width=600):
        def update(f=0):
            r.data_source.data['image'] = [v.imagemca[f]]
            r.data_source.data['x'] = [v.x_min[f]]
            r.data_source.data['y'] = [v.y_min[f]]
            r.data_source.data['dw'] = [v.x_max[f]-v.x_min[f]]
            r.data_source.data['dh'] = [v.y_max[f]-v.y_min[f]]

            push_notebook(handle=s)

        for i, val in enumerate(self.data):
            for k, v in val.items():
                p = figure(height=plot_height, width=plot_width, tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                           tools="pan,wheel_zoom,box_zoom,reset,hover,crosshair,save")
                p.x_range.range_padding = p.y_range.range_padding = 0

                # must give a vector of image data for image parameter
                color_mapper = LinearColorMapper(palette="Viridis256")

                simage = ColumnDataSource(data=dict(image=[v.imagemca[0]], x=[v.x_min[0]], y=[
                                          v.y_min[0]], dw=[v.x_max[0]-v.x_min[0]], dh=[v.y_max[0]-v.y_min[0]],))

                r = p.image(image='image', source=simage, x='x', y='y',
                            dw='dw', dh='dh', color_mapper=color_mapper, level="image")
                p.grid.grid_line_width = 0.5

                # Defining properties of color mapper
                color_bar = ColorBar(color_mapper=color_mapper,
                                     label_standoff=12,
                                     location=(0, 0),
                                     title='Counts')
                p.add_layout(color_bar, 'right')

                p.toolbar.logo = None

                if title != None:
                    p.title.text = str(title)
                else:
                    p.title.text = f'Image Movie for Scan {k}'
                if xlabel != None:
                    p.xaxis.axis_label = str(xlabel)
                else:
                    p.xaxis.axis_label = v.mcpRSXS_axes[0]
                if ylabel != None:
                    p.yaxis.axis_label = str(ylabel)
                else:
                    p.yaxis.axis_label = v.mcpRSXS_axes[1]

                s = show(p, notebook_handle=True)
                display(widgets.interact(update, f=(0, len(v.imagemca)-1)))

#########################################################################################


class LoadBeamline(Load1d):
    def load(self, basedir, file, key, **kwargs):
        """
        Load one or multiple specific scan(s) for selected streams.

        Parameters
        ----------
        basedir : string
            Specifiy the absolute or relative path to experimental data.
        file : string
            Specify the file name (either ASCII or HDF5).
        key : string
        **kwargs: multiple, optional
            Options:
                norm : boolean
                    Norm the spectra to [0,1].
                    default: True
                xoffset : list of tuples
                    Offset the x-axis by applying a polynomial fit.
                    default: None
                xcoffset : float
                    Offset x-axis by constant value.
                    default : None 
                yoffset : list of tuples
                    Offset the y-axis by applying a polynomial fit.
                    default : None 
                ycoffset : float
                    Offset y-axis by constant value.
                    default : None
        """

        # Append all REIXS scan objects to scan list in current object.
        self.data.append(loadSCAbeamline(basedir, file, key, **kwargs))
        self.type.append(key)
        self.x_stream.append('Scan Number')
        self.filename.append(file)

    def add(*args):
        raise UserWarning('Undefined')

    def subtract(*args):
        raise UserWarning('Undefined')


def getBL(basedir, file, stream, *args):
    get_single_beamline_value(basedir, file, stream, *args)
