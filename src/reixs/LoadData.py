# Scientific Modules
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, interp2d

from bokeh.io import push_notebook
from bokeh.plotting import show, figure
from bokeh.models import ColumnDataSource, HoverTool, LinearColorMapper, ColorBar, Span, Label

# Utilities
import os
from collections import defaultdict

# Widgets
import ipywidgets as widgets
from IPython.display import display
from ipyfilechooser import FileChooser

# Edge Dict
from .edges import EdgeDict

from .util import COLORP, all_list_entries_equal
from .xeol import *
from .add_subtract import ScanAddition, ScanSubtraction
from .sca import loadSCAscans
from .mca import loadMCAscans
from .mesh import loadMeshScans
from .rsxs_mcp import loadRSXS1dROIscans, loadRSXS2dROIscans, loadRSXSImageStack

#########################################################################################
#########################################################################################


class Load1d:
    def __init__(self):
        self.data = list()
        self.type = list()
        self.x_stream = list()
        self.filename = list()
        self.plot_lim_x = [":",":"]
        self.plot_lim_y = [":",":"]
        self.legend_loc = 'outside'
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()

    def load(self, basedir, file, x_stream, y_stream, *args, norm=True, is_XAS=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        self.data.append(loadSCAscans(basedir, file, x_stream, y_stream, *args, norm=norm,
                         is_XAS=is_XAS, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background, deriv=deriv, energyloss=energyloss, grid_x=grid_x))
        self.type.append(y_stream)
        self.x_stream.append(x_stream)
        self.filename.append(file)

    def add(self, basedir, file, x_stream, y_stream, *args, avg=False, norm=False, is_XAS=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        self.data.append(ScanAddition(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, is_XAS=is_XAS, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x))
        self.x_stream.append(x_stream)
        self.type.append(y_stream)
        self.filename.append(file)

    def subtract(self, basedir, file, x_stream, y_stream, *args, avg=False, norm=False, is_XAS=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        self.data.append(ScanSubtraction(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, is_XAS=is_XAS, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x))
        self.x_stream.append(x_stream)
        self.type.append(y_stream)
        self.filename.append(file)

    def xlim(self,lower,upper):
        self.plot_lim_x[0] = lower
        self.plot_lim_x[1] = upper

    def ylim(self,lower,upper):
        self.plot_lim_y[0] = lower
        self.plot_lim_y[1] = upper

    def plot_legend(self,pos):
        self.legend_loc = pos

    def vline(self,pos,**kwargs):
        self.plot_vlines.append([pos,kwargs])
    
    def hline(self,pos,**kwargs):
        self.plot_hlines.append([pos,kwargs])

    def label(self,pos_x,pos_y,text,**kwargs):
        self.plot_labels.append([pos_x,pos_y,text,kwargs])

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

        p = figure(plot_height=450,plot_width=700,
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

        if len(self.plot_hlines)>0:
            for line_props in self.plot_hlines:
                line = Span(location=line_props[0],dimension='width',**line_props[1])
                p.add_layout(line)

        if len(self.plot_vlines)>0:
            for line_props in self.plot_vlines:
                line = Span(location=line_props[0],dimension='height',**line_props[1])
                p.add_layout(line)

        if len(self.plot_labels)>0:
            for label_props in self.plot_labels:
                label = Label(x=label_props[0],y=label_props[1],text=label_props[2],**label_props[3])
                p.add_layout(label)

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
                df = df.append(s2)

        dfT = df.transpose(copy=True)

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


class XASLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=True, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None, deriv=None, grid_x=[None, None, None]):
        x_stream = 'Mono Energy'
        super().load(basedir, file, x_stream, y_stream, *args, norm=norm,
                     is_XAS=True, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background, deriv=deriv, grid_x=grid_x)

    def plot(self, linewidth=4):
        title = 'Absorption spectra normalized by mesh current'
        xlabel = "Incident Photon Energy (eV)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel)

    def add(self, basedir, file, y_stream, *args, avg=True, norm=True, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, grid_x=[None, None, None]):
        x_stream = "Mono Energy"
        super().add(basedir, file, x_stream, y_stream, *args, avg=avg, norm=norm,
                    is_XAS=True, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, grid_x=grid_x)

    def subtract(self, basedir, file, y_stream, *args, avg=True, norm=True, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, grid_x=[None, None, None]):
        x_stream = "Mono Energy"
        super().subtract(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, is_XAS=True, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, grid_x=grid_x)


class XESLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        x_stream = "MCP Energy"
        super().load(basedir, file, x_stream, y_stream, *
                     args, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x)

    def plot(self, linewidth=4):
        title = 'Summed MCP emission spectra'
        xlabel = "Uncalibrated MCP Energy (eV)"
        ylabel = "Counts (arb. units)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel, ylabel=ylabel)

    def add(self, basedir, file, y_stream, *args, avg=False, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        x_stream = "MCP Energy"
        super().add(basedir, file, x_stream, y_stream, *args,
                    avg=avg, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x)

    def subtract(self, basedir, file, y_stream, *args, avg=False, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        x_stream = "MCP Energy"
        super().subtract(basedir, file, x_stream, y_stream, *
                         args, avg=avg, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x)


class XRFLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        x_stream = "SDD Energy"
        super().load(basedir, file, x_stream, y_stream, *
                     args, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x)

    def plot(self, linewidth=4):
        title = 'Summed SDD emission spectra'
        xlabel = "Uncalibrated SDD Energy (eV)"
        ylabel = "Counts (arb. units)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel, ylabel=ylabel)

    def add(self, basedir, file, y_stream, *args, avg=False, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        x_stream = "SDD Energy"
        super().add(basedir, file, x_stream, y_stream, *args,
                    avg=avg, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x)

    def subtract(self, basedir, file, y_stream, *args, avg=False, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, energyloss=None, grid_x=[None, None, None]):
        x_stream = "SDD Energy"
        super().subtract(basedir, file, x_stream, y_stream, *
                         args, avg=avg, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, energyloss=energyloss, grid_x=grid_x)


class XEOLLoader(Load1d):
    def load(self, basedir, file, y_stream, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None, deriv=None, grid_x=[None, None, None]):
        x_stream = "XEOL Energy"
        super().load(basedir, file, x_stream, y_stream, *args, norm=norm,
                     xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background, deriv=deriv, grid_x=grid_x)

    def plot(self, linewidth=4):
        title = 'Summed XEOL spectra'
        xlabel = "Uncalibrated XEOL Wavelength (nm)"
        ylabel = "Counts (arb. units)"
        super().plot(linewidth=linewidth, title=title, xlabel=xlabel, ylabel=ylabel)

    def add(self, basedir, file, y_stream, *args, avg=False, norm=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, grid_x=[None, None, None]):
        x_stream = "XEOL Energy"
        super().add(basedir, file, x_stream, y_stream, *args, avg=avg,
                    norm=norm, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, grid_x=grid_x)

    def subtract(self, basedir, file, y_stream, *args, avg=False, norm=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None, grid_x=[None, None, None]):
        x_stream = "XEOL Energy"
        super().subtract(basedir, file, x_stream, y_stream, *args, avg=avg,
                         norm=norm, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv, grid_x=grid_x)


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
        self.grid_x = list()
        self.grid_y = list()
        self.plot_lim_x = [":",":"]
        self.plot_lim_y = [":",":"]
        self.plot_vlines = list()
        self.plot_hlines = list()
        self.plot_labels = list()


    def load(self, basedir, file, x_stream, y_stream, detector, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None, grid_x=[None, None, None], grid_y=[None, None, None]):
        if len(args) != 1:
            raise TypeError("You may only select one scan at a time")
        if self.data != []:
            raise TypeError("You can only append one scan per object")
        self.data.append(loadMCAscans(basedir, file, x_stream, y_stream, detector, *args, norm=norm,
                         xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background))
        self.x_stream.append(x_stream)
        self.y_stream.append(y_stream)
        self.detector.append(detector)
        self.norm.append(norm)
        self.filename.append(file)
        self.grid_x.append(grid_x)
        self.grid_y.append(grid_y)

    def xlim(self,lower,upper):
        self.plot_lim_x[0] = lower
        self.plot_lim_x[1] = upper

    def ylim(self,lower,upper):
        self.plot_lim_y[0] = lower
        self.plot_lim_y[1] = upper

    def vline(self,pos,**kwargs):
        self.plot_vlines.append([pos,kwargs])
    
    def hline(self,pos,**kwargs):
        self.plot_hlines.append([pos,kwargs])

    def label(self,pos_x,pos_y,text,**kwargs):
        self.plot_labels.append([pos_x,pos_y,text,kwargs])

    def plot(self):

        for i, val in enumerate(self.data):
            for k, v in val.items():

                p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                           tools="pan,wheel_zoom,box_zoom,reset,hover,crosshair,save")
                p.x_range.range_padding = p.y_range.range_padding = 0

                xmin, xmax, ymin, ymax, new_x, new_y, new_z = self.grid_data(
                    v, i)

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

                if self.plot_lim_y[0] != ':':
                    p.y_range.start = self.plot_lim_y[0]
                if self.plot_lim_y[1] != ':':
                    p.y_range.end = self.plot_lim_y[1]

                if self.plot_lim_x[0] != ':':
                    p.x_range.start = self.plot_lim_x[0]
                if self.plot_lim_x[1] != ':':
                    p.x_range.end = self.plot_lim_x[1]

                if len(self.plot_hlines)>0:
                    for line_props in self.plot_hlines:
                        line = Span(location=line_props[0],dimension='width',**line_props[1])
                        p.add_layout(line)

                if len(self.plot_vlines)>0:
                    for line_props in self.plot_vlines:
                        line = Span(location=line_props[0],dimension='height',**line_props[1])
                        p.add_layout(line)

                if len(self.plot_labels)>0:
                    for label_props in self.plot_labels:
                        label = Label(x=label_props[0],y=label_props[1],text=label_props[2],**label_props[3])
                        p.add_layout(label)

                p.title.text = f'{self.detector[i]} Image for Scan {k}'
                p.toolbar.logo = None
                p.xaxis.axis_label = str(self.x_stream[i])
                p.yaxis.axis_label = f"{self.y_stream[i]}"

                show(p)

    def grid_data(self, v, i):

        # Do auto-grid if not specified otherwise
        if self.grid_x[i] == [None, None, None]:
            xmin = v.x_data.min()
            xmax = v.x_data.max()
            x_points = int(
                np.ceil((xmax-xmin)/np.abs(np.diff(v.x_data)).min())) + 1

        else:
            xmin = self.grid_x[i][0]
            xmax = self.grid_x[i][1]
            x_points = int(np.ceil((xmax-xmin)/self.grid_x[i][2])) + 1

        if self.grid_y[i] == [None, None, None]:
            ymin = v.y_data.min()
            ymax = v.y_data.max()
            y_points = int(
                np.ceil((ymax-ymin)/np.abs(np.diff(v.y_data)).min())) + 1

        else:
            ymin = self.grid_y[i][0]
            ymax = self.grid_y[i][1]
            y_points = int(np.ceil((ymax-ymin)/self.grid_y[i][2])) + 1

        f = interp2d(v.x_data, v.y_data, v.detector)

        new_x = np.linspace(xmin, xmax, x_points, endpoint=True)
        new_y = np.linspace(ymin, ymax, y_points, endpoint=True)
        new_z = f(new_x, new_y)

        return xmin, xmax, ymin, ymax, new_x, new_y, new_z

    def export(self, filename):
        df = pd.DataFrame()
        with open(f"{filename}.txt_scale", "w") as f:
            with open(f"{filename}.txt_matrix", "w") as g:
                for i, val in enumerate(self.data):
                    for k, v in val.items():
                        xmin, xmax, ymin, ymax, new_x, new_y, new_z = self.grid_data(
                            v, i)
                        f.write("========================\n")
                        f.write(
                            f"F~{self.filename[i]}_S{v.scan}_{self.detector[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                        f.write("========================\n")

                        g.write("========================\n")
                        g.write(
                            f"F~{self.filename[i]}_S{v.scan}_{self.detector[i]}_{self.x_stream[i]}_{self.y_stream[i]}\n")
                        g.write("========================\n")

                        s1 = pd.Series(
                            new_x, name="Motor Scale Gridded")
                        df = df.append(s1)
                        s2 = pd.Series(
                            new_y, name="Detector Scale Gridded")
                        df = df.append(s2)
                        dfT = df.transpose(copy=True)
                        dfT.to_csv(f, index=False, line_terminator='\n')

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
    def __init__(self):
        self.energyloss = list()
        super().__init__()

    def load(self, basedir, file, detector, *args, norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None, grid_x=[None, None, None], grid_y=[None, None, None],energyloss=False):
        x_stream = 'Mono Energy'

        if detector == "MCP":
            y_stream = "MCP Energy"
        elif detector == "SDD":
            y_stream = "SDD Energy"
        elif detector == "XEOL":
            y_stream = "XEOL Energy"
        else:
            raise TypeError("Detector not defined.")

        self.energyloss.append(energyloss)

        super().load(basedir, file, x_stream, y_stream, detector, *args, norm=norm, xoffset=xoffset,
                     xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background, grid_x=grid_x, grid_y=grid_y)

    def grid_data(self, v, i):
        if self.energyloss[i] == False:
            return super().grid_data(v, i)
        
        # We will be swicthing x and y later
        # Here, grid_x and grid_y will already refer to switched axes
        else:
            if self.grid_y[i] == [None,None,None]:
                xmin = v.x_data.min()
                xmax = v.x_data.max()
                x_points = int(np.ceil((xmax-xmin)/np.abs(np.diff(v.x_data)).min())) + 1
            else:
                xmin = self.grid_y[i][0]
                xmax = self.grid_y[i][1]
                x_points = int(np.ceil((xmax-xmin)/self.grid_y[i][2])) + 1

            energy_loss_axes = list()

            for monoE in v.x_data:
                energy_loss_axes.append(monoE-v.y_data)

            if self.grid_x[i] == [None,None,None]:
                ymin = energy_loss_axes[-1][-1]
                ymax = energy_loss_axes[0][0]
                y_points = int(np.abs(np.ceil((ymax-ymin)/np.diff(energy_loss_axes[0]).min())))

            else:
                ymin = self.grid_x[i][0]
                ymax = self.grid_x[i][1]
                y_points = int(np.ceil((ymax-ymin)/self.grid_x[i][2])) + 1

            new_x = np.linspace(xmin, xmax, x_points, endpoint=True)
            new_y = np.linspace(ymin, ymax, y_points, endpoint=True)

            scatter_z = np.zeros((len(v.x_data),len(new_y)))

            for idx,val in enumerate(np.transpose(v.detector)):
                scatter_z[idx,:] = interp1d(energy_loss_axes[idx],val)(new_y)

            f = interp2d(v.x_data,new_y,np.transpose(scatter_z))
            new_z = f(new_x,new_y)

            # Switch x and y, so that we plot MCP Energy loss on horizontal axis
            # Overwrite stream names
            # Possibe since there can only be one scan per loader

            self.x_stream[i] = "Energy loss (eV)"
            self.y_stream[i] = "Mono Energy (eV)"

            return ymin, ymax, xmin, xmax, new_x, new_y, np.transpose(new_z)

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

#########################################################################################


class ImageROILoader():
    def __init__(self):
        self.data = list()
        self.filename = list()
        self.norm = list()

    def load(self, basedir, file, images, axis, *args, x=[None, None], y=[None, None], norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None):
        self.data.append(loadRSXS1dROIscans(basedir, file, images, axis, *args, x=x, y=y,
                         norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, deriv=deriv))
        self.norm.append(norm)
        self.filename.append(file)

    def plot(self, linewidth=4, title='Image ROI', xlabel=None, ylabel='Counts (arb. units)'):
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

        p = figure(plot_height=400,
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

    def export(self, filename):
        df = pd.DataFrame()
        files = list()

        for i, val in enumerate(self.data):
            for k, v in val.items():
                name = f"~{self.filename[i]}"
                if name not in files:
                    files.append(name)
                fileindex = files.index(name)

                for image, imagevalue in v.caxis.items():
                    s1 = pd.Series(
                        imagevalue, name=f"F{fileindex+1}_S{v.scan}_I{i+1}_Img{image}-{v.caxis_label}")
                    df = df.append(s1)
                    s2 = pd.Series(
                        v.imagesca[image], name=f"F{fileindex+1}_S{v.scan}_I{i+1}_Img{image}_Counts")
                    df = df.append(s2)

        dfT = df.transpose(copy=True)

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

    def load(self, basedir, file, axis, *args, x=[None, None], y=[None, None], norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None):
        self.data.append(loadRSXS2dROIscans(basedir, file, axis, *args, x=x, y=y, norm=norm,
                         xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset))
        self.norm.append(norm)
        self.filename.append(file)

    def plot(self):
        for i, val in enumerate(self.data):
            for k, v in val.items():

                p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
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

                p.title.text = f'Image Stack ROI for Scan {k}'
                p.toolbar.logo = None
                p.xaxis.axis_label = v.caxis_labels[0]
                p.yaxis.axis_label = v.caxis_labels[1]

                show(p)

    def export(self, filename):
        with open(f"{filename}.txt_scale", "a") as f:
            with open(f"{filename}.txt_matrix", "a") as g:
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

    def plot(self):
        def update(f=0):
            r.data_source.data['image'] = [v.imagemca[f]]
            r.data_source.data['x'] = [v.x_min[f]]
            r.data_source.data['y'] = [v.y_min[f]]
            r.data_source.data['dw'] = [v.x_max[f]-v.x_min[f]]
            r.data_source.data['dh'] = [v.y_max[f]-v.y_min[f]]

            push_notebook(handle=s)

        for i, val in enumerate(self.data):
            for k, v in val.items():
                p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
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

                p.title.text = f'Image Movie for Scan {k}'
                p.toolbar.logo = None
                p.xaxis.axis_label = v.mcpRSXS_axes[0]
                p.yaxis.axis_label = v.mcpRSXS_axes[1]

                s = show(p, notebook_handle=True)
                display(widgets.interact(update, f=(0, len(v.imagemca)-1)))
