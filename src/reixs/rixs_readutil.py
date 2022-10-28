# Imports
import numpy as np
import warnings

# Edge Dict
from .edges import EdgeDict

# Util functions
from .util import check_idx

def XES(mcp_data):
    """Sum the MCP detector image over all recorded datapoints."""
    try:
        sum = np.sum(mcp_data, axis=1)

    except:
        warnings.warn("Only one data point in MCP file to sum.")
        sum = mcp_data

    return sum


def rXES(mono_energy, mcp_data_norm, xes_incident_start, xes_incident_end):
    """Calculate resonant emission at selected energy (MCP)."""
    idx_start = (np.abs(xes_incident_start - mono_energy)).argmin()
    idx_end = (np.abs(xes_incident_end - mono_energy)).argmin()

    idx_start, idx_end = check_idx(idx_start,idx_end)

    return mcp_data_norm[idx_start:idx_end, :].sum(axis=0)


def XRF(sdd_data):
    """Sum the SDD detector image over all recorded datapoints"""
    return np.sum(sdd_data, axis=0)


def rXRF(mono_energy, sdd_data_norm, xes_incident_start, xes_incident_end):
    """Calculate resonant emission at selected energy (SDD)."""
    idx_start = (np.abs(xes_incident_start - mono_energy)).argmin()
    idx_end = (np.abs(xes_incident_end - mono_energy)).argmin()

    idx_start, idx_end = check_idx(idx_start,idx_end)

    return sdd_data_norm[idx_start:idx_end, :].sum(axis=0)


def ROI(PFY_edge=None, iPFY_edge=None):
    """Set a region of interest (ROI) over which to sum PFY (and iPFY)

        Parameters:
            PFY_edge : specify the absorption egde over which to sum.
                        Options:
                            N, O

            iPFY_edge : iPFY_edge when setting an ROI for iPFY
    """

    if PFY_edge != None:
        SDDLowerBound = EdgeDict[PFY_edge][0]
        SDDUpperBound = EdgeDict[PFY_edge][1]

        return SDDLowerBound, SDDUpperBound

    elif iPFY_edge != None:
        iSDDLowerBound = EdgeDict[iPFY_edge][0]
        iSDDUpperBound = EdgeDict[iPFY_edge][1]

        return iSDDLowerBound, iSDDUpperBound

    else:
        raise Exception("Proper ROI bounds not defined.")


def PFY(sdd_energy, sdd_data, mesh_current, PFY_edge=None, SDDLowerBound=None, SDDUpperBound=None):
    """ Calculate PFY based on a defined ROI.

    Prerequisites:
        Option 1: Define ROI via ROI method
        Option 2: Define ROI explictily by calling
                    SDDLowerBound : lower energy limit for summation
                    SDDUpperBound : upper energy limit for summation

    """
    if PFY_edge != None:
        SDDLowerBound, SDDUpperBound = ROI(PFY_edge)

    elif SDDLowerBound != None and SDDUpperBound != None:
        pass

    else:
        raise Exception(
            "Proper PFY boundaries not defined. Define boundaries or use TFY.")

    idx_low_sdd = (np.abs(SDDLowerBound - sdd_energy)).argmin()
    idx_high_sdd = (np.abs(SDDUpperBound - sdd_energy)).argmin()

    idx_low_sdd, idx_high_sdd = check_idx(idx_low_sdd,idx_high_sdd)

    # Transform this to PFY spectrum (rows = incident energy points, columns = detector 'space')
    sdd_detector_sum = sdd_data[:,
                                idx_low_sdd:idx_high_sdd].sum(axis=1)
    PFY_spec = (sdd_detector_sum/mesh_current)

    return PFY_spec


def iPFY(sdd_energy, sdd_data, mesh_current, iPFY_edge=None, iSDDLowerBound=None, iSDDUpperBound=None):
    """ Calculate iPFY based on a defined ROI.

        Prerequisites:
            Option 1: Define ROI via ROI method (and specify **kwarg iPFY_edge)
            Option 2: Define ROI explictily by calling
                    iSDDLowerBound : lower energy limit for summation
                    iSDDUpperBound : upper energy limit for summation

    """
    if iPFY_edge != None:
        iSDDLowerBound, iSDDUpperBound = ROI(iPFY_edge)

    elif iSDDLowerBound != None and iSDDUpperBound != None:
        pass

    else:
        raise Exception("Proper iPFY boundaries not defined.")

    iidx_low_sdd = (np.abs(iSDDLowerBound - sdd_energy)).argmin()
    iidx_high_sdd = (np.abs(iSDDUpperBound - sdd_energy)).argmin()

    iidx_low_sdd, iidx_high_sdd = check_idx(iidx_low_sdd,iidx_high_sdd)

    sdd_detector_sum_i = sdd_data[:,
                                  iidx_low_sdd:iidx_high_sdd].sum(axis=1)
    iPFY_spec = 1/(sdd_detector_sum_i/mesh_current)

    return iPFY_spec


def TFY(sdd_data, mesh_current):
    """Calculate TFY by summing over entire SDD image without ROI."""
    # Transform this to TFY spectrum (rows = incident energy points, columns = detector 'space')
    sdd_detector_sum = sdd_data[:, :].sum(axis=1)
    TFY_spec = (sdd_detector_sum/mesh_current)

    return TFY_spec


def specPFY(mcp_energy, mcp_data_norm, mcp_lowE, mcp_highE):
    """Calculate spectrometer PFY based on ROI set."""
    mcp_lowE_idx = (np.abs(mcp_lowE-mcp_energy)).argmin()
    mcp_highE_idx = (np.abs(mcp_highE-mcp_energy)).argmin()

    mcp_lowE_idx, mcp_highE_idx = check_idx(mcp_lowE_idx,mcp_highE_idx)

    return mcp_data_norm[:, mcp_lowE_idx:mcp_highE_idx].sum(axis=1)

def detector_norm(detector_mca,mesh_current):
    """Normalize the detector data by mesh current."""
    return np.true_divide(detector_mca,mesh_current)