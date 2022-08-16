import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def apply_offset(stream, offset=None, coffset=None):
    if offset != None:
        offsetarray = np.array(offset)
        coeff = np.polyfit(
            offsetarray[:, 0], offsetarray[:, 1], deg=len(offsetarray)-1)

        # Make sure that constant shift (only one tuple provided is handled as offset)
        if len(coeff) == 1:
            shift = offsetarray[0, 1] - offsetarray[0, 0]
            stream = stream+shift
        else:
            stream = np.polyval(coeff, stream)

    else:
        pass

    if coffset != None:
        return stream+coffset
    else:
        return stream

def grid_data(x_stream, y_stream, grid):
    xmin = grid[0]
    xmax = grid[1]

    numPoints = int(np.ceil((xmax-xmin)/grid[2])) + 1
    new_x = np.linspace(xmin, xmax, numPoints)
    f = interp1d(x_stream, y_stream, fill_value='extrapolate')
    new_y = f(new_x)

    return new_x, new_y

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

def apply_savgol(x,y,window,polyorder,deriv):
    xmin = x.min()
    xmax = x.max()
    x_diff = np.abs(np.diff(x)).min()
    new_x, new_y  = grid_data(x,y,[xmin,xmax,x_diff])

    if deriv == 0:
        delta = 1
    else:
        delta = x_diff
    smooth_y = savgol_filter(new_y,window,polyorder,deriv,delta)

    return new_x,smooth_y