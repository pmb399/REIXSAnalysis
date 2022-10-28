import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def apply_offset(stream, offset=None, coffset=None):
    """Apply constant or polynomial offset to specified stream
    
        Parameters
        ----------
        stream : array
            Specify the data to act on
        offset : list of tuples
            List all tuples with shift values (is,should)
        cofset : float
            Shift by constant
    """

    # Do the polynomial fitting with deg = len(list)-1
    if offset != None:
        offsetarray = np.array(offset)
        coeff = np.polyfit(
            offsetarray[:, 0], offsetarray[:, 1], deg=len(offsetarray)-1)

        # Make sure that constant shift is applied if only one tuple provided (handled as offset)
        if len(coeff) == 1:
            shift = offsetarray[0, 1] - offsetarray[0, 0]
            stream = stream+shift
        else:
            stream = np.polyval(coeff, stream)

    else:
        pass

    # Apply constant offset
    if coffset != None:
        return stream+coffset
    else:
        return stream

def grid_data(x_stream, y_stream, grid):
    """Grid 1d data
    
        Parameters
        ----------
        x_stream : array
            Specify the x data to act on 
        y_stream : array
            Specify the y data to act on
        grid : list, len 3
            Specify start value, end value, and delta
    """

    xmin = grid[0]
    xmax = grid[1]

    # Calculate number of data points
    numPoints = int(np.ceil((xmax-xmin)/grid[2])) + 1
    
    # Create linear space
    new_x = np.linspace(xmin, xmax, numPoints)

    # Do the interpolation step
    f = interp1d(x_stream, y_stream, fill_value='extrapolate')
    new_y = f(new_x)

    return new_x, new_y

def bin_data(x_data,y_data,binsize):
    """Reduce noise by averaging data points via binning mechanisms
    
        Parameters
        ----------
        x_data : array
            Specify the x data to act on
        y_data : array
            Specify the y data to act on
        binsize : int
            Specify how many data points to combine
            Must be exponents of 2
    """

    # Check if binsize is exponent of 2
    if (np.log(binsize)/np.log(2)).is_integer():
        pass
    else:
        raise ValueError("Bin size must be exponents of 2.")

    # Caluclate how many bins
    bins = len(x_data)/binsize

    # Split the data into the bins
    x_splits = np.split(x_data,bins)
    y_splits = np.split(y_data,bins)

    new_x = list()
    new_y = list()

    # Calculate the mean for all x and y values, respecitively, in the bin
    for idx,val in enumerate(x_splits):
        new_x.append(np.mean(val))
        new_y.append(np.mean(y_splits[idx]))

    return np.array(new_x), np.array(new_y)

### This is obsolete ###
def take_derivative1d(x, y, deg):
    if deg == 1:
        return np.gradient(y, x)
    elif deg == 2:
        dy = np.gradient(y, x)
        return np.gradient(dy, x)
    else:
        raise TypeError("No other derivatives implemented.")


def take_derivative2d(z, x, y, deg):
    if deg == 1:
        return np.gradient(z, x, y)
    elif deg == 2:
        dz = np.gradient(z, x, y)
        return np.gradient(dz, x, y)
    else:
        raise TypeError("No other derivatives implemented.")

### ---- ###

# New way for smooting and taking derivatives
def apply_savgol(x,y,window,polyorder,deriv):
    """Appply smoothing and take derivatives
    
        Parameters
        ----------
        x : array
            x data
        y : array
            y data
        window : int
            Length of the moving window in the Savitzgy-Golay filter
        polyorder : int
            Order of the fitted polynomial
        deriv : int
            Order of the derivative
            Choose "0" if only smoothing requested
    """
    xmin = x.min()
    xmax = x.max()
    # Caluclate minimum distance between data points to grid evenly spaced
    x_diff = np.abs(np.diff(x)).min()
    new_x, new_y  = grid_data(x,y,[xmin,xmax,x_diff])

    # Set default parameters as per savgol docs
    if deriv == 0:
        delta = 1
    else:
        delta = x_diff
    smooth_y = savgol_filter(new_y,window,polyorder,deriv,delta)

    return new_x,smooth_y