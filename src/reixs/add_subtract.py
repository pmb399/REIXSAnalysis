import numpy as np
from scipy.interpolate import interp1d
from .sca import loadSCAscans
from .simplemath import apply_offset, apply_savgol

def ScanAddition(basedir, file, x_stream, y_stream, *args, avg=True, norm=False, is_XAS=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None,energyloss=None,grid_x=[None,None,None],savgol=None,binsize=None):
    class added_object:
        def __init__(self):
            pass

    for i in args:
        if args.count(i) > 1:
            raise ValueError("Cannot add the same scan to itself")

    # Get the appropriate data first
    Scandata = loadSCAscans(basedir, file, x_stream, y_stream, *args,
                            norm=False, is_XAS=is_XAS, background=background,energyloss=None,grid_x=grid_x,binsize=binsize)

    for i, (k, v) in enumerate(Scandata.items()):
        if i == 0:
            MASTER_x_stream = v.x_stream
            MASTER_y_stream = v.y_stream
            name = str(k)+'+'
        else:
            if y_stream == 'XES' or y_stream.startswith('rXES'):
                if not np.array_equal(MASTER_x_stream, v.x_stream):
                    raise ValueError(
                        "Cannot add emission spectra with different energy scales.")
                else:
                    MASTER_y_stream += v.y_stream
            else:
                interp = interp1d(v.x_stream, v.y_stream,
                                  fill_value='extrapolate')(MASTER_x_stream)
                MASTER_y_stream += interp

            name += "_" + str(k)

    if avg == True:
        MASTER_y_stream = MASTER_y_stream/(i+1)

    data = dict()
    data[0] = added_object()
    data[0].x_stream = MASTER_x_stream
    data[0].y_stream = MASTER_y_stream
    data[0].scan = name

    if norm == True:
        data[0].y_stream = np.interp(
            data[0].y_stream, (data[0].y_stream.min(), data[0].y_stream.max()), (0, 1))

    data[0].x_stream = apply_offset(data[0].x_stream, xoffset, xcoffset)
    data[0].y_stream = apply_offset(data[0].y_stream, yoffset, ycoffset)

    if savgol != None:
        if isinstance(savgol,tuple):
            if len(savgol) == 2:
                savgol_deriv = 0
            elif len(savgol) == 3:
                savgol_deriv = savgol[2]
            else:
                raise TypeError("Savgol smoothing arguments incorrect.")
            data[0].x_stream, data[0].y_stream = apply_savgol(data[0].x_stream,data[0].y_stream,savgol[0],savgol[1],savgol_deriv)

            if norm == True:
                data[0].y_stream = data[0].y_stream / \
                data[0].y_stream.max()
        else:
            raise TypeError("Savgol smoothing arguments incorrect.")

    if energyloss!=None:
        data[0].x_stream = energyloss-data[0].x_stream

    return data


def ScanSubtraction(basedir, file, x_stream, y_stream, *args, avg=True, norm=False, is_XAS=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None,energyloss=None,grid_x=[None,None,None], savgol=None,binsize=None):
    class added_object:
        def __init__(self):
            pass

    for i in args:
        if args.count(i) > 1:
            raise ValueError("Cannot add the same scan to itself")

    # Get the appropriate data first
    Scandata = loadSCAscans(basedir, file, x_stream, y_stream, *args,
                            norm=False, is_XAS=is_XAS, background=background,energyloss=None,grid_x=grid_x,binsize=binsize)

    for i, (k, v) in enumerate(Scandata.items()):
        if i == 0:
            MASTER_x_stream = v.x_stream
            MASTER_y_stream = v.y_stream
            name = str(k) + '-'
        else:
            if y_stream == 'XES' or y_stream.startswith('rXES'):
                if not np.array_equal(MASTER_x_stream, v.x_stream):
                    raise ValueError(
                        "Cannot subtract emission spectra with different energy scales.")
                else:
                    MASTER_y_stream -= v.y_stream
            else:
                interp = interp1d(v.x_stream, v.y_stream,
                                  fill_value='extrapolate')(MASTER_x_stream)
                MASTER_y_stream -= interp

            name += "_" + str(k)

    if avg == True:
        MASTER_y_stream = MASTER_y_stream/(i+1)

    data = dict()
    data[0] = added_object()
    data[0].x_stream = MASTER_x_stream
    data[0].y_stream = MASTER_y_stream
    data[0].scan = name

    if norm == True:
        data[0].y_stream = np.interp(
            data[0].y_stream, (data[0].y_stream.min(), data[0].y_stream.max()), (0, 1))

    data[0].x_stream = apply_offset(data[0].x_stream, xoffset, xcoffset)
    data[0].y_stream = apply_offset(data[0].y_stream, yoffset, ycoffset)

    if savgol != None:
        if isinstance(savgol,tuple):
            if len(savgol) == 2:
                savgol_deriv = 0
            elif len(savgol) == 3:
                savgol_deriv = savgol[2]
            else:
                raise TypeError("Savgol smoothing arguments incorrect.")
            data[0].x_stream, data[0].y_stream = apply_savgol(data[0].x_stream,data[0].y_stream,savgol[0],savgol[1],savgol_deriv)

            if norm == True:
                data[0].y_stream = data[0].y_stream / \
                data[0].y_stream.max()
        else:
            raise TypeError("Savgol smoothing arguments incorrect.")

    if energyloss!=None:
        data[0].x_stream = energyloss-data[0].x_stream

    return data
