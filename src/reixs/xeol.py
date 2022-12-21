import numpy as np
from .util import check_idx

#########################################################################################


def XEOL_background_removal(data, arg, REIXSobj, background_scan=None):
    """Removes the background from XEOL scan
    
        Options:
            background_scan : True
                Uses the getBackground function in spec to subtract the background
            background_scan : int
                Uses a dedicated xeol scan in the current file to subtract
    """

    # Only apply background subtraction if specified
    if background_scan != None:
        # Determine where to load the background from
        if background_scan == True:
            background_data = data[arg].xeol_background
        elif type(background_scan) in [int]:
            background_data = REIXSobj.Scan(background_scan).xeol_data
        else:
            raise TypeError("Background Scan undefined.")
        
        ## If the background data has not been pre-procssed before or if only one frame was collected
        if background_data.ndim > 1:
            # 1 Sum the background
            background_sum = background_data.sum(axis=0)

            # 2 Normalize to average background frame
            background_spec = np.true_divide(background_sum, len(background_data))

        elif background_data.ndim == 1:
            background_spec = background_data

        else:
            raise TypeError("Background Scan has the wrong format")

        # 3 Subtract the background from the actual spectrum
        xeol_data = data[arg].xeol_data
        # To do so, need to broadcast background into shape
        xeol_subtracted = np.subtract(xeol_data, background_spec)

        return xeol_subtracted

    else:
        return data[arg].xeol_data


def xeol_spec(data, arg, REIXSobj, background_scan=None):
    """Calculate a XEOL spectrum on wavelength scale"""
    xeol_subtracted = XEOL_background_removal(
        data, arg, REIXSobj, background_scan)
    return np.sum(xeol_subtracted, axis=0)


def xeol_idx_lambda(data, arg, lambda_low, lambda_high):
    """Determine the indices associated wit lower and upper bounds on the xeol wavelength scale"""
    idx_low = (np.abs(lambda_low - data[arg].xeol_energy)).argmin()
    idx_high = (np.abs(lambda_high - data[arg].xeol_energy)).argmin()

    idx_low, idx_high = check_idx(idx_low,idx_high)

    return idx_low, idx_high


def xeol_idx_energy(data, arg, energy_low, energy_high):
    """Determine the indices assosciated with lower and upper bounds on the monochromator"""
    idx_low = (np.abs(energy_low - data[arg].mono_energy)).argmin()
    idx_high = (np.abs(energy_high - data[arg].mono_energy)).argmin()

    if idx_low == idx_high:
        idx_high = idx_low+1

    return idx_low, idx_high


def rxeol_spec(data, arg, REIXSobj, energy_low, energy_high, background_scan=None):
    """Calculates a XEOL spectrum as specific monochromator energy in an spec Escan/rscan"""
    # Remove background if requested
    xeol_subtracted = XEOL_background_removal(
        data, arg, REIXSobj, background_scan)
    # Get the indices of the monochromator window as specified
    ixd_low, idx_high = xeol_idx_energy(data, arg, energy_low, energy_high)
    # Ensure we sum and normalize by mesh
    return np.sum(np.true_divide(xeol_subtracted[ixd_low:idx_high, :], data[arg].mesh_current[ixd_low:idx_high, None]), axis=0)


def toy_spec(data, arg, REIXSobj, background_scan=None):
    """Calculates the total optical yield via integration over detector for every mono data point"""
    # Subtract background if requested
    xeol_subtracted = XEOL_background_removal(
        data, arg, REIXSobj, background_scan)
    # Ensure we sum and normalize by mesh
    return np.sum(np.true_divide(xeol_subtracted, data[arg].mesh_current[:, None]), axis=1)


def poy_spec(data, arg, REIXSobj, lambda_low, lambda_high, background_scan=None):
    """Calculates the partial optical yield given integration boundaries on wavelength scale"""
    # Subtract background if requested
    xeol_subtracted = XEOL_background_removal(
        data, arg, REIXSobj, background_scan)
    # Get indices of spectrometer window as specified by wavelengths
    ixd_low, idx_high = xeol_idx_lambda(data, arg, lambda_low, lambda_high)
    # Ensure we sum and normalize to flux
    return np.sum(np.true_divide(xeol_subtracted[:, ixd_low:idx_high], data[arg].mesh_current[:, None]), axis=1)


def xeol_map(data, arg, REIXSobj, background_scan=None):
    """Calculates a 2D XEOl map"""
    # Subtract background if requested
    xeol_subtracted = XEOL_background_removal(
        data, arg, REIXSobj, background_scan)
    # Normalize to incident flux (mesh)
    return np.true_divide(xeol_subtracted, data[arg].mesh_current[:, None])
