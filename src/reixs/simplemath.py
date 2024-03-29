import numpy as np
from scipy.interpolate import interp1d, interp2d
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
        # Limit the applicable shift to be quadratic or of lower order
        deg = min(2,len(offsetarray)-1)
        coeff = np.polyfit(
            offsetarray[:, 0], offsetarray[:, 1], deg=deg)

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

def grid_data2d(x_data, y_data, detector, grid_x=[None, None, None],grid_y=[None, None,None],energyloss=False):
    """Internal function to apply specified grid or ensure otherwise that axes are evenly spaced as this is required to plot an image."""
    if energyloss == False:
        # Do auto-grid if not specified otherwise
        # Take step-size as smallest delta observed in data array
        if grid_x == [None, None, None]:
            xmin = x_data.min()
            xmax = x_data.max()
            x_points = int(
                np.ceil((xmax-xmin)/np.abs(np.diff(x_data)).min())) + 1

        else:
            xmin = grid_x[0]
            xmax = grid_x[1]
            x_points = int(np.ceil((xmax-xmin)/grid_x[2])) + 1

        # Same as above, now for second axis.
        if grid_y == [None, None, None]:
            ymin = y_data.min()
            ymax = y_data.max()
            y_points = int(
                np.ceil((ymax-ymin)/np.abs(np.diff(y_data)).min())) + 1

        else:
            ymin = grid_y[0]
            ymax = grid_y[1]
            y_points = int(np.ceil((ymax-ymin)/grid_y[2])) + 1

        # Interpolate the data with given grid
        f = interp2d(x_data, y_data, detector)

        new_x = np.linspace(xmin, xmax, x_points, endpoint=True)
        new_y = np.linspace(ymin, ymax, y_points, endpoint=True)
        # Interpolate image on evenly-spaced grid
        new_z = f(new_x, new_y)
    
    elif energyloss == True:
        # We will be swicthing x and y later
        # Here, grid_x and grid_y will already refer to switched axes
        if grid_y == [None,None,None]:
            xmin = x_data.min()
            xmax = x_data.max()
            x_points = int(np.ceil((xmax-xmin)/np.abs(np.diff(x_data)).min())) + 1
        else:
            xmin = grid_y[0]
            xmax = grid_y[1]
            x_points = int(np.ceil((xmax-xmin)/grid_y[2])) + 1

        energy_loss_axes = list()

        # Calculate the energy loss axis for each mono energy as the incident energy will change the axis.
        for monoE in x_data:
            energy_loss_axes.append(monoE-y_data)

        # Determine the widest range where data is available on the rotated image.
        if grid_x == [None,None,None]:
            ymin = energy_loss_axes[-1][-1]
            ymax = energy_loss_axes[0][0]
            y_points = int(np.abs(np.ceil((ymax-ymin)/np.diff(energy_loss_axes[0]).min())))

        else:
            ymin = grid_x[0]
            ymax = grid_x[1]
            y_points = int(np.ceil((ymax-ymin)/grid_x[2])) + 1

        # Generate new axis as per the rotated image above.
        new_x = np.linspace(xmin, xmax, x_points, endpoint=True)
        new_y = np.linspace(ymin, ymax, y_points, endpoint=True)

        scatter_z = np.zeros((len(x_data),len(new_y)))

        # Evaluate the detector image on the new common energy axis
        for idx,val in enumerate(np.transpose(detector)):
            scatter_z[idx,:] = interp1d(energy_loss_axes[idx],val)(new_y)

        f = interp2d(x_data,new_y,np.transpose(scatter_z))
        new_z = f(new_x,new_y)

        # Switch x and y, so that we plot MCP Energy loss on horizontal axis
        # Overwrite stream names
        # Possibe since there can only be one scan per loader

        return ymin, ymax, xmin, xmax, new_x, new_y, np.transpose(new_z)
    else:
        raise ValueError('Could not determine if energyloss requested or not.')

    return xmin, xmax, ymin, ymax, new_x, new_y, new_z

def grid_data_mesh(x_data,y_data,z_data):
    """Internal function to generate scatter histogram for 3 independent SCA streams."""
    xmin = x_data.min()
    xmax = x_data.max()
    ymin = y_data.min()
    ymax = y_data.max()
    zmin = z_data.min()
    zmax = z_data.max()

    xunique = np.unique(x_data)
    yunique = np.unique(y_data)

    xbin = len(xunique)
    ybin = len(yunique)

    new_z, xedge, yedge = np.histogram2d(x_data, y_data, bins=[xbin, ybin], range=[
                                            [xmin, xmax], [ymin, ymax]], weights=z_data)
    new_z = np.transpose(new_z)

    return xmin, xmax, ymin, ymax, xedge, yedge, new_z, zmin, zmax

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

# New way for smoothing and taking derivatives
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