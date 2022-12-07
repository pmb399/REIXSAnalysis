import numpy as np
from scipy.interpolate import interp1d
from .sca import loadSCAscans
from .simplemath import apply_offset, apply_savgol
import warnings

def ScanAddition(basedir, file, x_stream, y_stream, *args, avg=True, norm=False, is_XAS=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None,energyloss=None,grid_x=[None,None,None],savgol=None,binsize=None,legend_item=None):
    """Internal function to handle scan addition.

        Parameters
        ----------
        args : Same as for the Load1d class
        kwargs: See Load1d class, but additional
        avg : boolean, optional
            Averages over the data after adding.
            default : True
    
    """

    # Define generic object in which all data will be stored
    class added_object:
        def __init__(self):
            pass

    # Ensure we only add a unique scan once
    for i in args:
        if args.count(i) > 1:
            raise ValueError("Cannot add the same scan to itself")

    # Get the appropriate data first - same loader as always
    Scandata = loadSCAscans(basedir, file, x_stream, y_stream, *args,
                            norm=False, is_XAS=is_XAS, background=background,energyloss=None,grid_x=grid_x,binsize=binsize)

    # Iterate over all loaded scans
    for i, (k, v) in enumerate(Scandata.items()):
        # Set the first scan as master data
        if i == 0:
            MASTER_x_stream = v.x_stream
            MASTER_y_stream = v.y_stream
            MASTER_mono = np.average(v.mono_energy)
            name = str(k)+'+'
        else:
            # Ensure that we only add emission scans when the spectrometer energy scale is identical
            if y_stream == 'XES' or y_stream.startswith('rXES'):
                if not np.array_equal(MASTER_x_stream, v.x_stream):
                    raise ValueError(
                        "Cannot add emission spectra with different energy scales.")
                else:
                    MASTER_y_stream += v.y_stream

                if MASTER_mono > 0.1+ np.average(v.mono_energy) or MASTER_mono < np.average(v.mono_energy) - 0.1:
                    warnings.warn("Changed average mono energy by more than 0.1 eV.")

            else:
                # For scans other than emission, set the first x-scale as master and interpolate all
                # data suczessively to ensure appropriate addition
                interp = interp1d(v.x_stream, v.y_stream,
                                  fill_value='extrapolate')(MASTER_x_stream)
                MASTER_y_stream += interp

            name += "_" + str(k)

    # Do the averaging here if requested.
    if avg == True:
        MASTER_y_stream = MASTER_y_stream/len(args)

    # Place data in a dictionary with the same structure as a regular Load1d call, so that we can plot it
    data = dict()
    data[0] = added_object()
    data[0].x_stream = MASTER_x_stream
    data[0].y_stream = MASTER_y_stream
    data[0].scan = name

    # Get legend items
    if legend_item != None:
        data[0].legend = legend_item
    else:
        data[0].legend = f"S{name}_{y_stream}"

    # Normalize data to [0,1]
    if norm == True:
        data[0].y_stream = np.interp(
            data[0].y_stream, (data[0].y_stream.min(), data[0].y_stream.max()), (0, 1))

    # May apply constant and polynomial offset
    data[0].x_stream = apply_offset(data[0].x_stream, xoffset, xcoffset)
    data[0].y_stream = apply_offset(data[0].y_stream, yoffset, ycoffset)

    # Apply smoothing and derivatives
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

    # Shift data to energy loss scale
    if energyloss!=None:
        if energyloss == True:
            data[0].x_stream = MASTER_mono-data[0].x_stream
        else:
            data[0].x_stream = energyloss-data[0].x_stream

    return data


def ScanSubtraction(basedir, file, x_stream, y_stream, *args, norm=False, is_XAS=False, background=None, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None,energyloss=None,grid_x=[None,None,None], savgol=None,binsize=None, legend_item=None):
    """Internal function to handle scan subtraction.
        New: May handle subtraction from two lists (addition within lists)

        Parameters
        ----------
        args : Same as for the Load1d class
            *args : ints comma separated or 2 lists
        kwargs: See Load1d class, but additional    
    """

    # Define generic object in which all data will be stored
    class added_object:
        def __init__(self):
            pass

    # Ensure we work with unique scans
    for i in args:
        if args.count(i) > 1:
            raise ValueError("Cannot subtract the same scan from itself")

    # Allows to define two lists
    # Add all scans within the two lists, then subtract results from each other
    if len(args) == 2 and type(args[0])==list and type(args[1])==list:
        minuend = ScanAddition(basedir, file, x_stream, y_stream, *args[0], avg=False, norm=norm, is_XAS=is_XAS, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset,energyloss=energyloss,grid_x=grid_x,savgol=savgol,binsize=binsize)
        subtrahend = ScanAddition(basedir, file, x_stream, y_stream, *args[1], avg=False, norm=norm, is_XAS=is_XAS, background=background, xoffset=xoffset, xcoffset=xcoffset, yoffset=yoffset, ycoffset=ycoffset,energyloss=energyloss,grid_x=grid_x,savgol=savgol,binsize=binsize)

        # Define the first scan (addition of scans as master)
        MASTER_x_stream = minuend[0].x_stream
        MASTER_y_stream = minuend[0].y_stream
        MASTER_mono = np.average(v.mono_energy)
        name = f"{args[0]}-{args[1]}"

        # Ensure same spectrometer energy scale between scans
        if y_stream == 'XES' or y_stream.startswith('rXES'):
            if not np.array_equal(MASTER_x_stream, subtrahend[0].x_stream):
                raise ValueError(
                    "Cannot subtract emission spectra with different energy scales.")
            else:
                MASTER_y_stream -= subtrahend[0].y_stream

            if MASTER_mono > 0.1+ np.average(v.mono_energy) or MASTER_mono < np.average(v.mono_energy) - 0.1:
                warnings.warn("Changed average mono energy by more than 0.1 eV.")

        else:
            # Interpolate all other data (not emission) to common x-scale.
            interp = interp1d(subtrahend[0].x_stream, subtrahend[0].y_stream,
                            fill_value='extrapolate')(MASTER_x_stream)
            MASTER_y_stream -= interp

    else:
        # Get the appropriate data first
        Scandata = loadSCAscans(basedir, file, x_stream, y_stream, *args,
                                norm=False, is_XAS=is_XAS, background=background,energyloss=None,grid_x=grid_x,binsize=binsize)

        # Iterate over all requested scans, load them, and subtract.
        # Same as above.
        for i, (k, v) in enumerate(Scandata.items()):
            if i == 0:
                MASTER_x_stream = v.x_stream
                MASTER_y_stream = v.y_stream
                MASTER_mono = np.average(v.mono_energy)
                name = str(k) + '-'
            else:
                if y_stream == 'XES' or y_stream.startswith('rXES'):
                    if not np.array_equal(MASTER_x_stream, v.x_stream):
                        raise ValueError(
                            "Cannot subtract emission spectra with different energy scales.")
                    else:
                        MASTER_y_stream -= v.y_stream

                    if MASTER_mono > 0.1+ np.average(v.mono_energy) or MASTER_mono < np.average(v.mono_energy) - 0.1:
                        warnings.warn("Changed average mono energy by more than 0.1 eV.")

                else:
                    interp = interp1d(v.x_stream, v.y_stream,
                                    fill_value='extrapolate')(MASTER_x_stream)
                    MASTER_y_stream -= interp

                name += "_" + str(k)

    # Place data in a dictionary with the same structure as a regular Load1d call, so that we can plot it
    data = dict()
    data[0] = added_object()
    data[0].x_stream = MASTER_x_stream
    data[0].y_stream = MASTER_y_stream
    data[0].scan = name

    # Get legend items
    if legend_item != None:
        data[0].legend = legend_item
    else:
        data[0].legend = f"S{name}_{y_stream}"

    # Normalize to [0,1]
    if norm == True:
        data[0].y_stream = np.interp(
            data[0].y_stream, (data[0].y_stream.min(), data[0].y_stream.max()), (0, 1))

    # Apply constant or polynomial offset.
    data[0].x_stream = apply_offset(data[0].x_stream, xoffset, xcoffset)
    data[0].y_stream = apply_offset(data[0].y_stream, yoffset, ycoffset)

    # Apply smoothing and derivative
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

    # Convert emission energy to energy loss scale if requested.
    if energyloss!=None:
        if energyloss == True:
            data[0].x_stream = MASTER_mono-data[0].x_stream
        else:
            data[0].x_stream = energyloss-data[0].x_stream

    return data
