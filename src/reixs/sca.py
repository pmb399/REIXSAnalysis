from .util import doesMatchPattern, check_key_in_dict
from .ReadData import REIXS
from .edges import EdgeDict
from .xeol import *
from .simplemath import apply_offset, take_derivative1d
import warnings
import numpy as np
from .parser import math_stream


def loadSCAscans(basedir, file, x_stream, y_stream, *args, norm=True, is_XAS=False, offset=None, coffset=None, background=None, deriv=None):
    special_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'XES', 'rXES', 'specPFY',
                       'XRF', 'rXRF', 'XEOL', 'rXEOL', 'POY', 'TOY', 'EY']  # all special inputs
    XAS_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'specPFY', 'POY',
                   'TOY', 'rXES', 'rXRF', 'rXEOL']  # All that are normalized to mesh

    # Note that the data dict only gets populated locally until the appropriate
    # singular y stream is return -- Assignment of y_stream happens after evaluation
    # of mathemtatical expression

    def get_y_data(y_stream, data, arg, background, REIXSObj):
        if doesMatchPattern(y_stream, special_streams):
            # CAUTION: Order matters as doesmatchpattern also returns TRUE if pattern in instance
            # is recogniced --> shortest pattern to check last

            if doesMatchPattern(y_stream, ['rXES']):
                roi = y_stream.lstrip("rXES[").rstrip("]")

                roi_low = int(roi.split(":")[0])
                roi_high = int(roi.split(":")[1])

                return data[arg].rXES(roi_low, roi_high)

            elif y_stream == 'XES':
                return data[arg].XES()

            elif doesMatchPattern(y_stream, ['rXRF']):
                roi = y_stream.lstrip("rXRF[").rstrip("]")

                roi_low = int(roi.split(":")[0])
                roi_high = int(roi.split(":")[1])
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
                    roi_low = int(roi.split(":")[0])
                    roi_high = int(roi.split(":")[1])
                    return data[arg].iPFY(iSDDLowerBound=roi_low, iSDDUpperBound=roi_high)

            elif doesMatchPattern(y_stream, ["specPFY"]):
                roi = y_stream.lstrip("specPFY[").rstrip("]")

                roi_low = int(roi.split(":")[0])
                roi_high = int(roi.split(":")[1])
                return data[arg].specPFY(roi_low, roi_high)

            elif doesMatchPattern(y_stream, ["PFY"]):
                roi = y_stream.lstrip("PFY[").rstrip("]")

                if check_key_in_dict(roi, EdgeDict):
                    return data[arg].PFY(PFY_edge=roi)

                else:
                    roi_low = int(roi.split(":")[0])
                    roi_high = int(roi.split(":")[1])
                    return data[arg].PFY(SDDLowerBound=roi_low, SDDUpperBound=roi_high)

            elif y_stream == 'XEOL':
                return xeol_spec(data, arg, REIXSobj, background_scan=background)

            elif doesMatchPattern(y_stream, ["rXEOL"]):
                roi = y_stream.lstrip("rXEOL[").rstrip("]")
                roi_low = int(roi.split(":")[0])
                roi_high = int(roi.split(":")[1])

                return rxeol_spec(data, arg, REIXSobj, roi_low, roi_high, background_scan=background)

            elif y_stream == 'TOY':
                return toy_spec(data, arg, REIXSobj, background_scan=background)

            elif doesMatchPattern(y_stream, ["POY"]):
                roi = y_stream.lstrip("POY[").rstrip("]")
                roi_low = int(roi.split(":")[0])
                roi_high = int(roi.split(":")[1])

                return poy_spec(data, arg, REIXSobj, roi_low, roi_high, background_scan=background)

            elif y_stream == 'EY':
                return data[arg].sample_current

            else:
                raise UserWarning("Special Stream not defined.")

        else:
            return np.array(data[arg].sca_data[y_stream])

    def get_x_data(x_stream, data, arg, background, REIXSObj):
        if x_stream == "Mono Energy":
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
            return np.array(data[arg].sca_data[x_stream])

    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        data[arg] = REIXSobj.Scan(arg)

        # Assign the calculated result to the y_stream of the object in data dict
        data[arg].y_stream = math_stream(
            y_stream, data, arg, get_y_data, XAS_streams=XAS_streams, is_XAS=is_XAS, background=background, REIXSObj=REIXSobj)

        # Assign the calculated result to the x_stream of the object in data dict
        data[arg].x_stream = math_stream(x_stream, data, arg, get_x_data)
        data[arg].x_stream = apply_offset(data[arg].x_stream, offset, coffset)

        if norm == True:
            data[arg].y_stream = np.interp(
                data[arg].y_stream, (data[arg].y_stream.min(), data[arg].y_stream.max()), (0, 1))

        if deriv != None:
            data[arg].y_stream = take_derivative1d(
                data[arg].y_stream, data[arg].x_stream, deriv)
            if norm == True:
                data[arg].y_stream = data[arg].y_stream / \
                    data[arg].y_stream.max()

    return data
