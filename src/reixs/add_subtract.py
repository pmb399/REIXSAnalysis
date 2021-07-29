import numpy as np
from scipy.interpolate import interp1d
from .sca import loadSCAscans
from .simplemath import apply_offset


def ScanAddition(basedir, file, x_stream, y_stream, *args, avg=True, norm=False, is_XAS=False, background=None, offset=None, coffset=None, deriv=None):
    class added_object:
        def __init__(self):
            pass

    for i in args:
        if args.count(i) > 1:
            raise ValueError("Cannot add the same scan to itself")

    # Get the appropriate data first
    Scandata = loadSCAscans(basedir, file, x_stream, y_stream, *args,
                            norm=False, is_XAS=is_XAS, background=background, deriv=deriv)

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

    data[0].x_stream = apply_offset(data[0].x_stream, offset, coffset)

    return data


def ScanSubtraction(basedir, file, x_stream, y_stream, *args, avg=True, norm=False, is_XAS=False, background=None, offset=None, coffset=None, deriv=None):
    class added_object:
        def __init__(self):
            pass

    for i in args:
        if args.count(i) > 1:
            raise ValueError("Cannot add the same scan to itself")

    # Get the appropriate data first
    Scandata = loadSCAscans(basedir, file, x_stream, y_stream, *args,
                            norm=False, is_XAS=is_XAS, background=background, deriv=deriv)

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

    data[0].x_stream = apply_offset(data[0].x_stream, offset, coffset)

    return data
