from .util import doesMatchPattern, check_key_in_dict, get_roi
from .ReadData import REIXS
from .edges import EdgeDict
from .xeol import *
from .simplemath import apply_offset, grid_data, apply_savgol, bin_data
import warnings
import numpy as np
from .parser import math_stream


def loadSCAscans(basedir, file, x_stream, y_stream, *args, norm=True, is_XAS=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None, energyloss=None, grid_x=[None, None, None], savgol=None, binsize=None):
    special_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'XES', 'rXES', 'specPFY',
                       'XRF', 'rXRF', 'XEOL', 'rXEOL', 'POY', 'TOY', 'EY', 'Sample', 'Mesh', 'ET']  # all special inputs
    XAS_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'specPFY', 'POY',
                   'TOY', 'rXES', 'rXRF', 'rXEOL', 'ET']  # All that are normalized to mesh

    # Note that the data dict only gets populated locally until the appropriate
    # singular y stream is return -- Assignment of y_stream happens after evaluation
    # of mathemtatical expression

    def get_y_data(y_stream, data, arg, background, REIXSObj):
        if doesMatchPattern(y_stream, special_streams):
            # CAUTION: Order matters as doesmatchpattern also returns TRUE if pattern in instance
            # is recogniced --> shortest pattern to check last

            if doesMatchPattern(y_stream, ['rXES']):
                roi = y_stream.lstrip("rXES[").rstrip("]")
                roi_low, roi_high = get_roi(roi)
                return data[arg].rXES(roi_low, roi_high)

            elif y_stream == 'XES':
                return data[arg].XES()

            elif doesMatchPattern(y_stream, ['rXRF']):
                roi = y_stream.lstrip("rXRF[").rstrip("]")
                roi_low, roi_high = get_roi(roi)
                return data[arg].rXRF(roi_low, roi_high)

            elif y_stream == 'XRF':
                return data[arg].XRF()

            elif y_stream == 'TEY':
                return data[arg].TEY

            elif y_stream == 'TFY':
                return data[arg].TFY()

            elif doesMatchPattern(y_stream, ["iPFY"]):
                roi = y_stream.lstrip("iPFY[").rstrip("]")

                if check_key_in_dict(roi, EdgeDict):
                    return data[arg].iPFY(iPFY_edge=roi)

                else:
                    roi_low, roi_high = get_roi(roi)
                    return data[arg].iPFY(iSDDLowerBound=roi_low, iSDDUpperBound=roi_high)

            elif doesMatchPattern(y_stream, ["specPFY"]):
                roi = y_stream.lstrip("specPFY[").rstrip("]")
                roi_low, roi_high = get_roi(roi)
                return data[arg].specPFY(roi_low, roi_high)

            elif doesMatchPattern(y_stream, ["PFY"]):
                roi = y_stream.lstrip("PFY[").rstrip("]")

                if check_key_in_dict(roi, EdgeDict):
                    return data[arg].PFY(PFY_edge=roi)

                else:
                    roi_low, roi_high = get_roi(roi)
                    return data[arg].PFY(SDDLowerBound=roi_low, SDDUpperBound=roi_high)

            elif y_stream == 'XEOL':
                return xeol_spec(data, arg, REIXSobj, background_scan=background)

            elif doesMatchPattern(y_stream, ["rXEOL"]):
                roi = y_stream.lstrip("rXEOL[").rstrip("]")
                roi_low, roi_high = get_roi(roi)

                return rxeol_spec(data, arg, REIXSobj, roi_low, roi_high, background_scan=background)

            elif y_stream == 'TOY':
                return toy_spec(data, arg, REIXSobj, background_scan=background)

            elif doesMatchPattern(y_stream, ["POY"]):
                roi = y_stream.lstrip("POY[").rstrip("]")
                roi_low, roi_high = get_roi(roi)

                return poy_spec(data, arg, REIXSobj, roi_low, roi_high, background_scan=background)

            if doesMatchPattern(y_stream, ['ET']):
                roi = y_stream.lstrip("ET[").rstrip("]")
                roi_low, roi_high = get_roi(roi)

                from .LoadData import EEMsLoader
                import pandas as pd
                import io

                eems = EEMsLoader()
                eems.load(basedir, file, 'MCP', arg, norm=norm, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset, background=background, grid_x=grid_x,energyloss=True)
                f,g = eems.get_data()

                df = pd.read_csv(io.StringIO(f.getvalue()),skiprows=3,delimiter=",")
                x_data = np.array(df["Motor Scale Gridded"].dropna())
                matrix = np.loadtxt(io.StringIO(g.getvalue()),skiprows=4)
                energy_transfer = np.array(df['Detector Scale Gridded'].dropna())
                
                idx_min = np.abs(roi_low-energy_transfer).argmin()
                idx_max = np.abs(roi_high-energy_transfer).argmin()

                y_data = np.sum(matrix[:,idx_min:idx_max],axis=1)
                data[arg].mono_energy_gridded = x_data

                return y_data

            elif y_stream == 'EY':
                return data[arg].sample_current

            # elif y_stream == 'tey':
            #    # spec mnemonic for sample current
            #    return data[arg].sample_current

            elif y_stream == 'Sample':
                # to ensure backwards compatibility
                return data[arg].sample_current

            elif y_stream == 'Mesh':
                # to ensure backwards compatibility
                return data[arg].mesh_current

            else:
                try:
                    return get_sca_data(y_stream, data, arg)

                except:
                    raise UserWarning("Special Stream not defined.")

        else:
            return get_sca_data(y_stream, data, arg)

    def get_sca_data(stream, data, arg):
        try:
            return np.array(data[arg].sca_data[stream])
        except:
            return np.array(data[arg].sca_data[data[arg].mnemonic2name[stream]])

    def get_x_data(x_stream, data, arg, background, REIXSObj):
        if x_stream == "Mono Energy":
            if y_stream.startswith("ET"):
                return data[arg].mono_energy_gridded
            else:
                return data[arg].mono_energy

        elif x_stream == "MCP Energy":
            return data[arg].mcp_energy

        elif x_stream == "SDD Energy":
            return data[arg].sdd_energy

        elif x_stream == "XEOL Energy":
            return data[arg].xeol_energy

        elif x_stream == "Points":
            return np.array(range(0, len(data[arg].y_stream)), dtype=int)

        else:
            return get_sca_data(x_stream, data, arg)

    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        data[arg] = REIXSobj.Scan(arg)

        # Assign the calculated result to the y_stream of the object in data dict
        data[arg].y_stream = math_stream(
            y_stream, data, arg, get_y_data, XAS_streams=XAS_streams, is_XAS=is_XAS, background=background, REIXSObj=REIXSobj)

        # Assign the calculated result to the x_stream of the object in data dict
        data[arg].x_stream = math_stream(x_stream, data, arg, get_x_data)

        #Bin the data if requested
        if binsize != None:
            data[arg].x_stream, data[arg].y_stream = bin_data(data[arg].x_stream,data[arg].y_stream,binsize)

        # Grid the data if specified
        if grid_x != [None, None, None]:
            new_x, new_y = grid_data(
                data[arg].x_stream, data[arg].y_stream, grid_x)

            data[arg].x_stream = new_x
            data[arg].y_stream = new_y

        # Apply offsets and normalize
        data[arg].x_stream = apply_offset(
            data[arg].x_stream, xoffset, xcoffset)

        if norm == True:
            data[arg].y_stream = np.interp(
                data[arg].y_stream, (data[arg].y_stream.min(), data[arg].y_stream.max()), (0, 1))

        data[arg].y_stream = apply_offset(
            data[arg].y_stream, yoffset, ycoffset)
               
        if savgol != None:
            if isinstance(savgol,tuple):
                if len(savgol) == 2:
                    savgol_deriv = 0
                elif len(savgol) == 3:
                    savgol_deriv = savgol[2]
                else:
                    raise TypeError("Savgol smoothing arguments incorrect.")
                data[arg].x_stream, data[arg].y_stream = apply_savgol(data[arg].x_stream,data[arg].y_stream,savgol[0],savgol[1],savgol_deriv)

                if norm == True:
                    data[arg].y_stream = data[arg].y_stream / \
                    data[arg].y_stream.max()
            else:
                raise TypeError("Savgol smoothing arguments incorrect.")

        # Transforms RIXS to energy loss scale if incident energy is given
        if energyloss != None:
            if energyloss == True:
                data[arg].x_stream = np.average(data[arg].mono_energy)-data[arg].x_stream
            else:
                data[arg].x_stream = energyloss-data[arg].x_stream

    return data
