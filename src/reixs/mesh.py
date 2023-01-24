from .util import doesMatchPattern, check_key_in_dict, get_roi
from .ReadData import REIXS
from .edges import EdgeDict
from .xeol import *
from .simplemath import apply_offset, grid_data_mesh
import warnings
import numpy as np
from .parser import math_stream


def loadMeshScans(basedir, file, x_stream, y_stream, z_stream, *args, norm=True, is_XAS=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None):
    """Internal function to generate scatter plots for (x,y,z) SCA data"""

    # Define special streams one might want to consider
    special_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'XES', 'rXES', 'specPFY',
                       'XRF', 'rXRF', 'XEOL', 'rXEOL', 'POY', 'TOY', 'EY']  # all special inputs
    # Append also to XAS_Streams variable, if the special stream is already normalized to i0 (incident flux)
    XAS_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'specPFY', 'POY',
                   'TOY', 'rXES', 'rXRF', 'rXEOL']  # All that are normalized to mesh

    # Note that the data dict only gets populated locally until the appropriate
    # singular y stream is return -- Assignment of y_stream happens after evaluation
    # of mathemtatical expression

    def get_z_data(z_stream, data, arg, background, REIXSObj):
        # Get the z-data stream
        # Check if the requested stream is a special stream
        if doesMatchPattern(z_stream, special_streams):
            # CAUTION: Order matters as doesmatchpattern also returns TRUE if pattern in instance
            # is recogniced --> shortest pattern to check last

            # Check which special stream is requested and strip the arguments for ROI if applicable
            if doesMatchPattern(z_stream, ['rXES']):
                roi = z_stream.lstrip("rXES[").rstrip("]")

                roi_low, roi_high = get_roi(roi)

                return data[arg].rXES(roi_low, roi_high)

            elif z_stream == 'XES':
                return data[arg].XES()

            elif doesMatchPattern(z_stream, ['rXRF']):
                roi = z_stream.lstrip("rXRF[").rstrip("]")
                roi_low, roi_high = get_roi(roi)

                return data[arg].rXRF(roi_low, roi_high)

            elif z_stream == 'XRF':
                return data[arg].XRF()

            elif z_stream == 'TEY':
                return data[arg].TEY

            elif z_stream == 'TFY':
                return data[arg].TFY()

            elif doesMatchPattern(z_stream, ["iPFY"]):
                roi = z_stream.lstrip("iPFY[").rstrip("]")

                if check_key_in_dict(roi, EdgeDict):
                    return data[arg].iPFY(iPFY_edge=roi)

                else:
                    roi_low, roi_high = get_roi(roi)
                    return data[arg].iPFY(iSDDLowerBound=roi_low, iSDDUpperBound=roi_high)

            elif doesMatchPattern(z_stream, ["specPFY"]):
                roi = z_stream.lstrip("specPFY[").rstrip("]")
                roi_low, roi_high = get_roi(roi)
                return data[arg].specPFY(roi_low, roi_high)

            elif doesMatchPattern(z_stream, ["PFY"]):
                roi = z_stream.lstrip("PFY[").rstrip("]")

                if check_key_in_dict(roi, EdgeDict):
                    return data[arg].PFY(PFY_edge=roi)

                else:
                    roi_low, roi_high = get_roi(roi)
                    return data[arg].PFY(SDDLowerBound=roi_low, SDDUpperBound=roi_high)

            elif z_stream == 'XEOL':
                return xeol_spec(data, arg, REIXSobj, background_scan=background)

            elif doesMatchPattern(z_stream, ["rXEOL"]):
                roi = z_stream.lstrip("rXEOL[").rstrip("]")
                roi_low, roi_high = get_roi(roi)

                return rxeol_spec(data, arg, REIXSobj, roi_low, roi_high, background_scan=background)

            elif z_stream == 'TOY':
                return toy_spec(data, arg, REIXSobj, background_scan=background)

            elif doesMatchPattern(z_stream, ["POY"]):
                roi = z_stream.lstrip("POY[").rstrip("]")
                roi_low, roi_high = get_roi(roi)

                return poy_spec(data, arg, REIXSobj, roi_low, roi_high, background_scan=background)

            elif z_stream == 'EY':
                return data[arg].sample_current

            else:
                try:
                    # Else, load from pandas SCA data frame
                    return get_sca_data(z_stream, data, arg)
                except:
                    raise UserWarning("Special Stream not defined.")

        else:
            return get_sca_data(z_stream, data, arg)

    def get_sca_data(stream, data, arg):
        # Read in stream with long name or utilize mnemonic dict
        try:
            return np.array(data[arg].sca_data[stream])
        except:
            raise UserWarning("Stream not defined. Only mnemonics supported!")

    def get_x_data(x_stream, data, arg, background, REIXSObj):
        # Get the x-data.
        if x_stream == "Mono Energy":
            return data[arg].mono_energy

        elif x_stream == "MCP Energy":
            return data[arg].mcp_energy

        elif x_stream == "SDD Energy":
            return data[arg].sdd_energy

        elif x_stream == "XEOL Energy":
            return data[arg].xeol_energy

        else:
            return get_sca_data(x_stream, data, arg)

    def get_y_data(y_stream, data, arg, background, REIXSObj):
        # Get the y data
        if y_stream == "Mono Energy":
            return data[arg].mono_energy

        elif y_stream == "MCP Energy":
            return data[arg].mcp_energy

        elif y_stream == "SDD Energy":
            return data[arg].sdd_energy

        elif y_stream == "XEOL Energy":
            return data[arg].xeol_energy

        else:
            return get_sca_data(y_stream, data, arg)

    # Generate dictionary to store data with REIXS objects
    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        # Load scans to dict
        data[arg] = REIXSobj.Scan(arg)

        # Assign the calculated result to the y_stream of the object in data dict
        data[arg].y_data = math_stream(y_stream, data, arg, get_y_data)
        data[arg].y_data = apply_offset(data[arg].y_data, yoffset, ycoffset)

        # Assign the calculated result to the x_stream of the object in data dict
        data[arg].x_data = math_stream(x_stream, data, arg, get_x_data)
        data[arg].x_data = apply_offset(data[arg].x_data, xoffset, xcoffset)

        # Aplly simple math to x-stream
        data[arg].z_data = math_stream(z_stream, data, arg, get_z_data, XAS_streams=XAS_streams,
                                       is_XAS=is_XAS, background=background, REIXSObj=REIXSobj)

        # Normalize if requested
        if norm == True:
            data[arg].z_data = np.interp(
                data[arg].z_data, (data[arg].z_data.min(), data[arg].z_data.max()), (0, 1))

        xmin, xmax, ymin, ymax, xedge, yedge, new_z, zmin, zmax = grid_data_mesh(data[arg].x_data,data[arg].y_data,data[arg].z_data)
        data[arg].xmin = xmin
        data[arg].xmax = xmax
        data[arg].ymin = ymin
        data[arg].ymax = ymax
        data[arg].xedge = xedge
        data[arg].yedge = yedge
        data[arg].new_z = new_z
        data[arg].zmin = zmin
        data[arg].zmax = zmax

    return data
