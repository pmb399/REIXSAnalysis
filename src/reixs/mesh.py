from reixs.util import doesMatchPattern, check_key_in_dict
from reixs.ReadData import REIXS
from reixs.edges import EdgeDict
from reixs.xeol import *
from reixs.offset import apply_offset
import warnings
import numpy as np
import re
import parser


def loadMeshScans(basedir, file, x_stream, y_stream, z_stream, *args, norm=True, is_XAS=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, background=None):
    special_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'XES', 'rXES', 'specPFY',
                       'XRF', 'rXRF', 'XEOL', 'rXEOL', 'POY', 'TOY']  # all special inputs
    XAS_streams = ['TEY', 'TFY', 'PFY', 'iPFY', 'specPFY', 'POY',
                   'TOY', 'rXES', 'rXRF', 'rXEOL']  # All that are normalized to mesh

    # Note that the data dict only gets populated locally until the appropriate
    # singular y stream is return -- Assignment of y_stream happens after evaluation
    # of mathemtatical expression

    def get_z_data(y_stream, data, arg, background, REIXSObj):
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

        else:
            return np.array(data[arg].sca_data[x_stream])

    def get_y_data(y_stream, data, arg, background, REIXSObj):
        if y_stream == "Mono Energy":
            return data[arg].mono_energy

        elif y_stream == "MCP Energy":
            return data[arg].mcp_energy

        elif y_stream == "SDD Energy":
            return data[arg].sdd_energy

        elif y_stream == "XEOL Energy":
            return data[arg].xeol_energy

        else:
            return np.array(data[arg].sca_data[y_stream])


    def math_stream(formula, data, arg, get_data, XAS_streams=None, is_XAS=False, background=None, REIXSObj=None):
        # Split the user input string at all mathematical operations
        # Allow "( ) * / + -" as math

        pattern = '[\(+\-*^/\)]'
        split_expr = re.split(pattern, formula)

        # Place string literals in dict if cannot be converted to float
        # drop all empty strings from re.split
        quantity_str_dict = dict()

        for i, string in enumerate(split_expr):
            if string != "":
                try:
                    float(string)
                except:
                    # Assign generic "val{i}" key to string literal in compliance with
                    # python supported syntax for variables
                    quantity_str_dict[f"val{i}"] = string

        # Parser does not support special string literals (ROIs) due to inproper python variable naming syntax
        # Replace them with generic key as per dictionary
        # Create local variable and assign corresponding data -- needed as eval interprets all input as variable
        # This is indeed local to this function only and cannot be accessed from loadSCAscans
        for k, v in quantity_str_dict.items():
            formula = formula.replace(v, k)

            # Ensure that in XASLoader all quantities are normalized by mesh
            # per definition of what XAS is
            # Exclude those quantities from normalization that have been normalized elsewhere

            if is_XAS == False:
                locals()[k] = get_data(v, data, arg, background, REIXSObj)
            else:
                if doesMatchPattern(v, XAS_streams):
                    locals()[k] = get_data(v, data, arg, background, REIXSObj)
                else:
                    numerator = get_data(v, data, arg, background, REIXSObj)
                    mesh = data[arg].mesh_current
                    locals()[k] = numerator/mesh

        # Return the calculated result
        code = parser.expr(formula).compile()
        return eval(code)

    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        data[arg] = REIXSobj.Scan(arg)

        # Assign the calculated result to the y_stream of the object in data dict
        data[arg].y_data = math_stream(y_stream, data, arg, get_y_data)
        data[arg].y_data = apply_offset(data[arg].y_data, yoffset, ycoffset)

        # Assign the calculated result to the x_stream of the object in data dict
        data[arg].x_data = math_stream(x_stream, data, arg, get_x_data)
        data[arg].x_data = apply_offset(data[arg].x_data, xoffset, xcoffset)

        data[arg].z_data = math_stream(z_stream, data, arg, get_z_data, XAS_streams=XAS_streams,
                                       is_XAS=is_XAS, background=background, REIXSObj=REIXSobj)

        if norm == True:
            data[arg].z_data = np.interp(
                data[arg].z_data, (data[arg].z_data.min(), data[arg].z_data.max()), (0, 1))

    return data
