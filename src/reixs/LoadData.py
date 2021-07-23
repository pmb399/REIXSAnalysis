# Scientific Modules
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, interp2d

from bokeh.io import show, output_notebook
from bokeh.plotting import show, figure
from bokeh.models import ColumnDataSource, HoverTool, LinearColorMapper, ColorBar
from bokeh.models.tools import LassoSelectTool, CrosshairTool

# Utilities
import os
import re
import parser
from collections import defaultdict
import warnings

# Widgets
import ipywidgets as widgets
from IPython.display import display, clear_output
from ipyfilechooser import FileChooser

# Edge Dict
from reixs.edges import EdgeDict

from reixs.ReadData import REIXS
from reixs.util import COLORP
from reixs.xeol import *
from reixs.math import ScanAddition, ScanSubtraction
from reixs.sca import loadSCAscans
from reixs.mca import loadMCAscans
from reixs.mesh import loadMeshScans

#########################################################################################
#########################################################################################


class Load1d:
    def __init__(self):
        self.data = list()
        self.type = list()
        self.x_stream = list()
        self.filename = list()

    def load(self, basedir, file, x_stream, y_stream, *args, norm=True, is_XAS=False, offset=None, coffset=None, background=None):
        self.data.append(loadSCAscans(basedir, file, x_stream, y_stream, *args, norm=norm,
                         is_XAS=is_XAS, offset=offset, coffset=coffset, background=background))
        self.type.append(y_stream)
        self.x_stream.append(x_stream)
        self.filename.append(file)

    def add(self, basedir, file, x_stream, y_stream, *args, avg=False, norm=False, is_XAS=False, background=None, offset=None, coffset=None):
        self.data.append(ScanAddition(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, is_XAS=is_XAS, background=background, offset=offset, coffset=coffset))
        self.x_stream.append(x_stream)
        self.type.append(y_stream)
        self.filename.append(file)

    def subtract(self, basedir, file, x_stream, y_stream, *args, avg=False, norm=False, is_XAS=False, background=None, offset=None, coffset=None):
        self.data.append(ScanSubtraction(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, is_XAS=is_XAS, background=background, offset=offset, coffset=coffset))
        self.x_stream.append(x_stream)
        self.type.append(y_stream)
        self.filename.append(file)

    def plot(self, linewidth=4, title=None, xlabel=None, ylabel=None):
        plot_data = defaultdict(list)
        for i, val in enumerate(self.data):
            for k, v in val.items():
                plot_data["x_stream"].append(v.x_stream)
                plot_data["y_stream"].append(v.y_stream)
                plot_data['x_name'].append(self.x_stream[i])
                plot_data['y_name'].append(self.type[i])
                plot_data['filename'].append(self.filename[i])
                plot_data['scan'].append(v.scan)
                plot_data['legend'].append(f"S{v.scan}_{self.type[i]}")

        numlines = len(plot_data['scan'])
        plot_data['color'] = COLORP[0:numlines]

        source = ColumnDataSource(plot_data)

        p = figure(plot_height=400,
                   tools="pan,wheel_zoom,box_zoom,reset,crosshair,save")
        p.multi_line(xs='x_stream', ys='y_stream', legend_group="legend",
                     line_width=linewidth, line_color='color', line_alpha=0.6,
                     hover_line_color='color', hover_line_alpha=1.0,
                     source=source)

        p.add_tools(HoverTool(show_arrow=False, line_policy='next', tooltips=[
            ('Scan', '@scan'),
            ('File', '@filename'),
            ("(x,y)", "(@x_name, @y_name)"),
            ("(x,y)", "($x, $y)")
        ]))

        p.toolbar.logo = None

        if title != None:
            p.title.text = str(title)
        if xlabel != None:
            p.xaxis.axis_label = str(xlabel)
        if ylabel != None:
            p.yaxis.axis_label = str(ylabel)
        show(p)

    def export(self, filename):
        df = pd.DataFrame()
        files = list()

        for i, val in enumerate(self.data):
            for k, v in val.items():
                name = f"~{self.filename[i]}"
                if name not in files:
                    files.append(name)
                fileindex = files.index(name)

                s1 = pd.Series(
                    v.x_stream, name=f"F{fileindex+1}_S{v.scan}_I{i+1}-{self.x_stream[i]}")
                df = df.append(s1)
                s2 = pd.Series(
                    v.y_stream, name=f"F{fileindex+1}_S{v.scan}_I{i+1}-{self.type[i]}")
                df = df.append(s2,)

        dfT = df.transpose(copy=True)

        with open(f"{filename}.csv", 'w') as f:
            for idx, file in enumerate(files):
                f.write(f"# F{idx+1} {file}\n")
            dfT.to_csv(f, index=False)

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


class XASLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=True, offset=None, coffset=None, background=None):
        x_stream = 'Mono Energy'
        super().load(basedir, file, x_stream, y_stream, *args, norm=norm,
                     is_XAS=True, offset=offset, coffset=coffset, background=background)

    def plot(self, linewidth=4):
        title = 'Absorption spectra normalized by mesh current'
        xlabel = "Incident Photon Energy (eV)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel)

    def add(self, basedir, file, y_stream, *args, avg=True, norm=True, background=None, offset=None, coffset=None):
        x_stream = "Mono Energy"
        super().add(basedir, file, x_stream, y_stream, *args, avg=avg, norm=norm,
                         is_XAS=True, background=background, offset=offset, coffset=coffset)

    def subtract(self, basedir, file, y_stream, *args, avg=True, norm=True, background=None, offset=None, coffset=None):
        x_stream = "Mono Energy"
        super().subtract(basedir, file, x_stream, y_stream, *args, avg=avg,
                              norm=norm, is_XAS=True, background=background, offset=offset, coffset=coffset)


class XESLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=False, offset=None, coffset=None):
        x_stream = "MCP Energy"
        super().load(basedir, file, x_stream, y_stream, *
                     args, norm=norm, offset=offset, coffset=coffset)

    def plot(self, linewidth=4):
        title = 'Summed MCP emission spectra'
        xlabel = "Uncalibrated MCP Energy (eV)"
        ylabel = "Counts (arb. units)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel, ylabel=ylabel)

    def add(self, basedir, file, y_stream, *args, avg=False, norm=False, offset=None, coffset=None):
        x_stream = "MCP Energy"
        super().add(basedir, file, x_stream, y_stream, *args,
                         avg=avg, norm=norm, offset=offset, coffset=coffset)

    def subtract(self, basedir, file, y_stream, *args, avg=False, norm=False, offset=None, coffset=None):
        x_stream = "MCP Energy"
        super().subtract(basedir, file, x_stream, y_stream, *
                              args, avg=avg, norm=norm, offset=offset, coffset=coffset)


class XRFLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=False, offset=None, coffset=None):
        x_stream = "SDD Energy"
        super().load(basedir, file, x_stream, y_stream, *
                     args, norm=norm, offset=offset, coffset=coffset)

    def plot(self, linewidth=4):
        title = 'Summed SDD emission spectra'
        xlabel = "Uncalibrated SDD Energy (eV)"
        ylabel = "Counts (arb. units)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel, ylabel=ylabel)

    def add(self, basedir, file, y_stream, *args, avg=False, norm=False, offset=None, coffset=None):
        x_stream = "SDD Energy"
        super().add(basedir, file, x_stream, y_stream, *args,
                         avg=avg, norm=norm, offset=offset, coffset=coffset)

    def subtract(self, basedir, file, y_stream, *args, avg=False, norm=False, offset=None, coffset=None):
        x_stream = "SDD Energy"
        super().subtract(basedir, file, x_stream, y_stream, *
                              args, avg=avg, norm=norm, offset=offset, coffset=coffset)


class XEOLLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=False, offset=None, coffset=None, background=None):
        x_stream = "XEOL Energy"
        super().load(basedir, file, x_stream, y_stream, *args, norm=norm,
                     offset=offset, coffset=coffset, background=background)

    def plot(self, linewidth=4):
        title = 'Summed XEOL spectra'
        xlabel = "Uncalibrated XEOL Wavelength (nm)"
        ylabel = "Counts (arb. units)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel, ylabel=ylabel)

    def add(self, basedir, file, y_stream, *args, avg=False, norm=False, background=None, offset=None, coffset=None):
        x_stream = "XEOL Energy"
        super().add(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, background=background, offset=offset, coffset=coffset)

    def subtract(self, basedir, file, y_stream, *args, avg=False, norm=False, background=None, offset=None, coffset=None):
        x_stream = "XEOL Energy"
        super().subtract(basedir, file, x_stream, y_stream, *args, avg=avg,
                              norm=norm, background=background, offset=offset, coffset=coffset)


#########################################################################################

class Load2d:
    """Plots the 2D image of a given detector"""

    def __init__(self):
        self.data = list()
        self.x_stream = list()
        self.y_stream = list()
        self.detector = list()
        self.filename = list()
        self.norm = list()

    def load(self, basedir, file, x_stream, y_stream, detector, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None):
        self.data.append(loadMCAscans(basedir, file, x_stream, y_stream, detector, *args, norm=norm,
                         xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background))
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.detector.append(detector)
        self.norm.append(norm)
        self.filename.append(file)

    def plot(self):

        for i, val in enumerate(self.data):
            for k, v in val.items():

                p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                           tools="pan,wheel_zoom,box_zoom,reset,hover,crosshair,save")
                p.x_range.range_padding = p.y_range.range_padding = 0

                xmin, xmax, ymin, ymax, new_x, new_y, new_z = self.grid_data(v)

                # must give a vector of image data for image parameter
                color_mapper = LinearColorMapper(palette="Viridis256",
                                                 low=v.detector.min(),
                                                 high=v.detector.max())

                p.image(image=[new_z], x=xmin, y=ymin, dw=xmax-xmin,
                        dh=ymax-ymin, color_mapper=color_mapper, level="image")
                p.grid.grid_line_width = 0.5

                # Defining properties of color mapper
                color_bar = ColorBar(color_mapper=color_mapper,
                                     label_standoff=12,
                                     location=(0, 0),
                                     title='Counts')
                p.add_layout(color_bar, 'right')

                if self.y_stream[i] == 'SDD Energy':
                    p.y_range.start = 0
                    p.y_range.end = 1200

                p.title.text = f'{self.detector[i]} Image for Scan {k}'
                p.toolbar.logo = None
                p.xaxis.axis_label = str(self.x_stream[i])
                p.yaxis.axis_label = f"{self.y_stream[i]}"

                show(p)

    def grid_data(self, v):
        xmin = v.x_data.min()
        xmax = v.x_data.max()
        ymin = v.y_data.min()
        ymax = v.y_data.max()

        x_points = int(np.ceil((xmax-xmin)/np.diff(v.x_data).min()))
        y_points = int(np.ceil((ymax-ymin)/np.diff(v.y_data).min()))

        f = interp2d(v.x_data, v.y_data, v.detector)

        new_x = np.linspace(xmin, xmax, x_points)
        new_y = np.linspace(ymin, ymax, y_points)
        new_z = f(new_x, new_y)

        return xmin, xmax, ymin, ymax, new_x, new_y, new_z

    def export(self, filename):
        with open(f"{filename}.txt_scale", "a") as f:
            with open(f"{filename}.txt_matrix", "a") as g:
                for i, val in enumerate(self.data):
                    for k, v in val.items():
                        xmin, xmax, ymin, ymax, new_x, new_y, new_z = self.grid_data(
                            v)
                        f.write("========================\n")
                        f.write(
                            f"F~{self.filename[i]}_S{v.scan}_{self.detector[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                        f.write("========================\n")

                        g.write("========================\n")
                        g.write(
                            f"F~{self.filename[i]}_S{v.scan}_{self.detector[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                        g.write("========================\n")

                        f.write("=== Motor Scale Gridded ===\n")
                        np.savetxt(f, new_x)
                        f.write("=== Detector Scale Gridded ===\n")
                        np.savetxt(f, new_y)
                        g.write("=== Image ===\n")
                        np.savetxt(g, new_z, fmt="%.9g")

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


class EEMsLoader(Load2d):
    def load(self, basedir, file, detector, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None):
        x_stream = 'Mono Energy'

        if detector == "MCP":
            y_stream = "MCP Energy"
        elif detector == "SDD":
            y_stream = "SDD Energy"
        elif detector == "XEOL":
            y_stream = "XEOL Energy"
        else:
            raise TypeError("Detector not defined.")

        super().load(basedir, file, x_stream, y_stream, detector, *args, norm=norm, xoffset=xoffset,
                     xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background)

#########################################################################################


class LoadMesh:
    """Plots the Mesh Scan images"""

    def __init__(self):
        self.data = list()
        self.x_stream = list()
        self.y_stream = list()
        self.z_stream = list()
        self.filename = list()
        self.norm = list()

    def load(self, basedir, file, x_stream, y_stream, z_stream, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None):
        self.data.append(loadMeshScans(basedir, file, x_stream, y_stream, z_stream, *args, norm=norm,
                         xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background))
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.z_stream.append(z_stream)
        self.norm.append(norm)
        self.filename.append(file)

    def plot(self):

        for i, val in enumerate(self.data):
            for k, v in val.items():

                p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                           tools="pan,wheel_zoom,box_zoom,reset,hover,crosshair,save")
                p.x_range.range_padding = p.y_range.range_padding = 0

                xmin, xmax, ymin, ymax, new_x, new_y, new_z, zmin, zmax = self.grid_data(
                    v)

                # must give a vector of image data for image parameter
                color_mapper = LinearColorMapper(palette="Viridis256",
                                                 low=zmin,
                                                 high=zmax)

                p.image(image=[new_z], x=xmin, y=ymin, dw=xmax-xmin,
                        dh=ymax-ymin, color_mapper=color_mapper, level="image")
                p.grid.grid_line_width = 0.5

                # Defining properties of color mapper
                color_bar = ColorBar(color_mapper=color_mapper,
                                     label_standoff=12,
                                     location=(0, 0),
                                     title='Counts')
                p.add_layout(color_bar, 'right')

                p.title.text = f'{self.z_stream[i]} Histogram for Scan {k}'
                p.toolbar.logo = None
                p.xaxis.axis_label = str(self.x_stream[i])
                p.yaxis.axis_label = f"{self.y_stream[i]}"

                show(p)

    def grid_data(self, v):
        xmin = v.x_data.min()
        xmax = v.x_data.max()
        ymin = v.y_data.min()
        ymax = v.y_data.max()
        zmin = v.z_data.min()
        zmax = v.z_data.max()

        xunique = np.unique(v.x_data)
        yunique = np.unique(v.y_data)

        xbin = len(xunique)
        ybin = len(yunique)

        new_z, xedge, yedge = np.histogram2d(v.x_data, v.y_data, bins=[xbin, ybin], range=[
                                             [xmin, xmax], [ymin, ymax]], weights=v.z_data)

        return xmin, xmax, ymin, ymax, xedge, yedge, new_z, zmin, zmax

    def export(self, filename):
        with open(f"{filename}.txt_scale", "a") as f:
            with open(f"{filename}.txt_matrix", "a") as g:
                for i, val in enumerate(self.data):
                    for k, v in val.items():
                        xmin, xmax, ymin, ymax, new_x, new_y, new_z, zmin, zmax = self.grid_data(
                            v)
                        f.write("========================\n")
                        f.write(
                            f"F~{self.filename[i]}_S{v.scan}_{self.z_stream[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                        f.write("========================\n")

                        g.write("========================\n")
                        g.write(
                            f"F~{self.filename[i]}_S{v.scan}_{self.z_stream[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                        g.write("========================\n")

                        f.write("=== x-axis bins ===\n")
                        np.savetxt(f, new_x)
                        f.write("=== y-axis bins ===\n")
                        np.savetxt(f, new_y)
                        g.write("=== Histogram ===\n")
                        np.savetxt(g, new_z, fmt="%.9g")

        print(f"Successfully wrote Histogram data to {filename}.txt")

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


class MeshLoader(LoadMesh):
    pass
