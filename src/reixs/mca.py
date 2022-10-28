import numpy as np
from .ReadData import REIXS
from .xeol import *
from .simplemath import apply_offset
from .parser import math_stream


def loadMCAscans(basedir, file, x_stream, y_stream, detector, *args, norm=True, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None):
    """Internal function to load MCA data
    
        Parameters
        ----------
        See Load2d function.
    """
    # Note that the data dict only gets populated locally until the appropriate
    # singular y stream is return -- Assignment of y_stream happens after evaluation
    # of mathemtatical expression

    def get_y_data(detector, data, arg, background, REIXSobj):
        # Sets the detector energy
        if detector == "SDD Energy":
            y_data = data[arg].sdd_energy
        elif detector == 'MCP Energy':
            y_data = data[arg].mcp_energy
        elif detector == "XEOL Energy":
            y_data = data[arg].xeol_energy

        else:
            raise TypeError("Stream undefined!")

        return y_data

    def get_z_data(detector, data, arg, background, REIXSobj):
        # Sets the detector MCA
        if detector == "SDD":
            z_data = np.transpose(data[arg].SDD_norm())
        elif detector == 'MCP':
            z_data = np.transpose(data[arg].MCP_norm())
        elif detector == 'XEOL':
            if background != None:
                z_data = np.transpose(
                    xeol_map(data, arg, REIXSobj, background_scan=background))
            else:
                z_data = np.transpose(xeol_map(data, arg, REIXSobj))

        else:
            raise TypeError("Detector undefined!")

        return z_data

    def get_x_data(x_stream, data, arg, background, REIXSobj):
        # Sets the independent axis (x-axis values)
        if x_stream == "Mono Energy":
            return data[arg].mono_energy

        else:
            try:
                return np.array(data[arg].sca_data[x_stream])
            except:
                return np.array(data[arg].sca_data[data[arg].mnemonic2name[x_stream]])

    # Place all loaded REIXS objects in data dictionary
    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        # Load all scans and place in data dict
        data[arg] = REIXSobj.Scan(arg)

        # Assign the calculated result to the y_stream of the object in data dict
        # May apply math operations and offsets
        data[arg].y_data = math_stream(y_stream, data, arg, get_y_data)
        data[arg].y_data = apply_offset(data[arg].y_data, yoffset, ycoffset)

        # Assign the calculated result to the x_stream of the object in data dict
        # May apply math operations and offsets
        data[arg].x_data = math_stream(x_stream, data, arg, get_x_data)
        data[arg].x_data = apply_offset(data[arg].x_data, xoffset, xcoffset)

        # Apply simple math operations on 2D MCA detector data
        data[arg].detector = math_stream(
            detector, data, arg, get_z_data, background=background, REIXSObj=REIXSobj)

        # Normalize if requested
        if norm == True:
            data[arg].detector = data[arg].detector/np.max(data[arg].detector)

    return data
