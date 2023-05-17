# Scientific Modules
import numpy as np
import pandas as pd
import h5py

# Edge Dict
from .edges import EdgeDict

# Scan analysis utils
from .rsxs_readutil import img_to_sca, grid_stack, stack_to_mca
from . import rixs_readutil

# Spec Config
from .spec_config import get_REIXSconfig

# Utilities
import os
import warnings

## Every scan will be an instance of a class ##

#########################################################################################


def REIXS(baseName, header_file):
    """Function that will load an REIXS data object either from ASCII or HDF5, 
    depending on extension of the header file

        Parameters
        ----------
        baseName : string
            Give directory path of stored data file
        header_file : string
            Give name of header file with extension
    """

    def handle_exception(e):
        return REIXS_ASCII(baseName, header_file, REIXSconfig)

    # Read in parameter/variable config (either external from os.env or internal default)
    REIXSconfig = get_REIXSconfig()

    # Return REIXS objects based on extension
    # Always default to HDF5
    if "." in header_file:
        return REIXS_HDF5(baseName, header_file, REIXSconfig, handle_exception)
    else:
        return REIXS_HDF5(baseName, str(header_file)+".h5", REIXSconfig, handle_exception)


class REIXS_HDF5(object):
    """REIXS data object  HDF5

        Parameters
        ----------
        baseName : string
            Give directory path of stored data file
        header_file : string
            Give name of header file with extension
        REIXSconfig : dict
            Spec paramter to python variable dict
    """

    def __init__(self, baseName, header_file, REIXSconfig, handle_exception):
        try:
            self.file = os.path.join(baseName, header_file)
        except:
            raise TypeError("You did not specify a directory path.")
        self.REIXSconfig = REIXSconfig

        self.handle_exception = handle_exception

    def Scan(self, scan):
        """Load one specific scan from specified data file

        Parameters
        ----------
        scan : int
        """
        return self.myScan(self, scan)

    class myScan(object):

        def __init__(my, self, scan):
            # Try opening hdf5 container
            try:
                with h5py.File(self.file, 'r') as f:

                    try:
                        # Create dictionary for scan numbers
                        scanIndexDictHeader = dict()
                        for k in f.keys():
                            # Groups must start with the SCAN prefix
                            if k.startswith("SCAN_"):
                                scanIndexDictHeader[int(k.split("_")[1])] = k
                        my.scan = scanIndexDictHeader[scan]

                    except:
                        raise KeyError("Scan not defined.")

                    # Define selected special streams and detectors
                    try:
                        my.mono_energy = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_mono_energy"]}'])
                    except:
                        raise UserWarning("Problem detecting energy.")

                    try:
                        my.mesh_current = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_mesh_current"]}'])
                    except:
                        raise UserWarning("Problem detecting mesh current.")

                    try:
                        my.sample_current = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_sample_current"]}'])
                    except:
                        raise UserWarning("Problem detecting sample current.")

                    try:
                        my.TEY = my.sample_current/my.mesh_current
                    except:
                        raise ValueError("Problem calculating TEY.")

                    try:
                        my.sdd_data = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_sdd_data"]}'])
                        my.sdd_energy = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_sdd_energy"]}'])
                    except:
                        warnings.warn(
                            "Could not load SDD data / SDD energy scale")

                    try:
                        my.xeol_data = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_xeol_data"]}'])
                        my.xeol_energy = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_xeol_energy"]}'])
                        my.xeol_background = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_xeol_background"]}'])

                    except:
                        warnings.warn(
                            "Could not load XEOL data / XEOL emission scale")

                    try:
                        my.mcp_data = np.transpose(
                            np.array(f[f'{my.scan}/{self.REIXSconfig["HDF5_mcp_data"]}']))
                        my.mcp_energy = np.array(
                            f[f'{my.scan}/{self.REIXSconfig["HDF5_mcp_energy"]}'])
                    except:
                        warnings.warn(
                            "Could not load MCP data / MCP energy scale")

                    # Populate a pandas dataframe with all SCA data
                    my.sca_data = pd.DataFrame()
                    try:
                        for entry in f[f'{my.scan}/{self.REIXSconfig["HDF5_sca_data"]}']:
                            if len(f[f'{my.scan}/{self.REIXSconfig["HDF5_sca_data"]}/{entry}'].shape) == 1 and len(f[f'{my.scan}/{self.REIXSconfig["HDF5_sca_data"]}/{entry}']) == len(f[f'{my.scan}/{self.REIXSconfig["HDF5_sca_data"]}/epoch']):
                                my.sca_data[str(entry)] = np.array(
                                    f[f'{my.scan}/{self.REIXSconfig["HDF5_sca_data"]}/{entry}'])
                    except:
                        warnings.warn(
                            "Could not load SCAs from HDF5 container.")

            except OSError as e:
                obj = self.handle_exception(e)

                my.scan = scan
                my.mono_energy = obj.Scan(scan).mono_energy
                my.mesh_current = obj.Scan(scan).mesh_current
                my.sample_current = obj.Scan(scan).sample_current
                my.TEY = my.sample_current/my.mesh_current
                try:
                    my.sdd_data = obj.Scan(scan).sdd_data
                    my.sdd_energy = obj.Scan(scan).sdd_energy
                except:
                    warnings.warn("Could not load SDD data / SDD energy scale")
                try:
                    my.xeol_data = obj.Scan(scan).xeol_data
                    my.xeol_energy = obj.Scan(scan).xeol_energy
                    my.xeol_background = obj.Scan(scan).xeol_background
                except:
                    warnings.warn(
                        "Could not load XEOL data / XEOL emission scale")
                try:
                    if obj.MCPRIXS == True:
                        my.mcp_data = obj.Scan(scan).mcp_data
                        my.mcp_energy = obj.Scan(scan).mcp_energy
                    elif obj.MCPRSXS == True:
                        my.mcpRSXS_scales = obj.Scan(scan).mcpRSXS_scales
                        my.mcpRSXS_scatters = obj.Scan(scan).mcpRSXS_scatters
                        my.mcpRSXS_axes = obj.Scan(scan).mcpRSXS_axes
                    else:
                        warnings.warn("No MCP data found")
                except:
                    warnings.warn("Could not load MCP data / MCP energy scale")
                try:
                    my.sca_data = obj.Scan(scan).sca_data
                except:
                    warnings.warn("Could not load SCAs from ASCII.")

            except:
                raise ValueError("Scan Data not defined")

        def MCP_norm(my):
            """Normalize the counts of the MCP by incident flux at every given datapoint.
               This is only applied to EEMs."""

            my.mcp_data_norm = np.transpose(
                rixs_readutil.detector_norm(my.mcp_data, my.mesh_current))

            return my.mcp_data_norm

        def SDD_norm(my):
            """Normalize the counts of the SDD by incident flux at every given datapoint."""
            my.sdd_data_norm = rixs_readutil.detector_norm(
                my.sdd_data, my.mesh_current[:, None])

            return my.sdd_data_norm

        def XES(my):
            """Sum the MCP detector image over all recorded datapoints."""

            return rixs_readutil.XES(my.mcp_data)

        def rXES(my, xes_incident_start, xes_incident_end):
            """Calculate resonant emission at selected energy (MCP)."""

            if not(hasattr(my, 'mcp_data_norm')):
                my.MCP_norm()

            return rixs_readutil.rXES(my.mono_energy, my.mcp_data_norm, xes_incident_start, xes_incident_end)

        def XRF(my):
            """Sum the SDD detector image over all recorded datapoints"""

            return rixs_readutil.XRF(my.sdd_data)

        def rXRF(my, xes_incident_start, xes_incident_end):
            """Calculate resonant emission at selected energy (SDD)."""

            if not(hasattr(my, 'sdd_data_norm')):
                my.SDD_norm()

            return rixs_readutil.rXRF(my.mono_energy, my.sdd_data_norm, xes_incident_start, xes_incident_end)

        def PFY(my, PFY_edge=None, SDDLowerBound=None, SDDUpperBound=None):
            """ Calculate PFY based on a defined ROI.

            Prerequisites:
                Option 1: Define ROI via ROI method
                    PFY_edge : string, name of element
                Option 2: Define ROI explictily by calling
                    SDDLowerBound : float, lower energy limit for summation
                    SDDUpperBound : float, upper energy limit for summation

            """

            return rixs_readutil.PFY(my.sdd_energy, my.sdd_data, my.mesh_current, PFY_edge, SDDLowerBound, SDDUpperBound)

        def iPFY(my, iPFY_edge=None, iSDDLowerBound=None, iSDDUpperBound=None):
            """ Calculate iPFY based on a defined ROI.

                Prerequisites:
                    Option 1: Define ROI via ROI method (and specify **kwarg iPFY_edge)
                        iPFY_edge : string, name of element
                    Option 2: Define ROI explictily by calling
                        iSDDLowerBound : float, lower energy limit for summation
                        iSDDUpperBound : float, upper energy limit for summation
            """

            return rixs_readutil.iPFY(my.sdd_energy, my.sdd_data, my.mesh_current, iPFY_edge, iSDDLowerBound, iSDDUpperBound)

        def TFY(my):
            """Calculate TFY by summing over entire SDD image without ROI."""
            # Transform this to TFY spectrum (rows = incident energy points, columns = detector 'space')

            return rixs_readutil.TFY(my.sdd_data, my.mesh_current)

        def specPFY(my, mcp_lowE, mcp_highE):
            """Calculate spectrometer PFY based on ROI set.

            Parameters
            ----------
            mcp_lowE : float
                Lower energy cutoff for detector integration
            mcp_highE : float
                Higher energy cutoff for detector integration
            """

            # Normalize MCP data to incident flux first, if not already done previously
            if not(hasattr(my, 'mcp_data_norm')):
                my.MCP_norm()

            return rixs_readutil.specPFY(my.mcp_energy, my.mcp_data_norm, mcp_lowE, mcp_highE)

        def RSXS_MCPnorm(my):
            """Normalize RSXS MCP data by incident flux"""
            my.RSXSMCP_norm = dict()
            # Need to appyly the flux correction to all scatters at each data point
            for k, v in my.mcpRSXS_scatters.items():
                my.RSXSMCP_norm[k] = np.true_divide(v, my.mesh_current[k])

            return my.RSXSMCP_norm

        def RSXS_1dROI(my, img, x_low=None, x_high=None, y_low=None, y_high=None, axis=0):
            """Collape RSXS MCP data to 1d plot"""
            if not(hasattr(my, 'RSXSMCP_norm')):
                my.RSXS_MCPnorm()

            return img_to_sca(my.mcpRSXS_scales, my.RSXSMCP_norm, img, x_low, x_high, y_low, y_high, axis)

        def RSXS_2dROI(my, x_low=None, x_high=None, y_low=None, y_high=None, axis=0):
            """Integrate RSXS MCP data to 2d plot"""
            if not(hasattr(my, 'RSXSMCP_norm')):
                my.RSXS_MCPnorm()

            return stack_to_mca(my.mcpRSXS_scales, my.RSXSMCP_norm, x_low, x_high, y_low, y_high, axis)

        def RSXS_Images(my):
            """Stack the data with multiple images."""
            if not(hasattr(my, 'RSXSMCP_norm')):
                my.RSXS_MCPnorm()

            return grid_stack(my.mcpRSXS_scales, my.RSXSMCP_norm)

    def Info(self, keys):
        """Load one specific scan from specified data file

        Parameters
        ----------
        keys : list of key paths to the desired information
        """
        # Try opening hdf5 container
        try:
            with h5py.File(self.file, 'r') as f:
                # Create dictionary for scan numbers
                info_dict = dict()

                if not isinstance(keys, list):
                    keys = [keys]

                for key in keys:
                    info_dict[key] = dict()
                    for k in f.keys():
                        if k.startswith("SCAN_"):  # Groups must start with the SCAN prefix
                            try:
                                info_dict[key][int(
                                    k.split("_")[1])] = f[f'{k}/{key}'][()].decode("utf-8")
                            except AttributeError:
                                entry = f[f'{k}/{key}'][()]
                                if isinstance(entry, np.ndarray) and len(entry)==1:
                                    entry = entry[0]
                                info_dict[key][int(
                                    k.split("_")[1])] = entry

        except:
            raise KeyError("Error opening and processing file.")

        return info_dict


class REIXS_ASCII(object):

    """Returns all scans in a scan file."""

    def __init__(self, baseName, header_file, REIXSconfig):
        """
        Constructor for one data file.

        Parameters
        ----------
        baseName : string
            Absolute path to experimental data.
        header_file : string
            Name of the header file.

        Returns
        -------
        REIXS : all data associated with the specified data file ready to be processed.
        """

        # Load all scans in specified file
        self.datasets = [[]]
        self.scanNumbers = []
        self.scanType = []
        self.scanMotor = []
        self.mnemonics_motors = []
        self.full_names_motors = []
        self.mnemonics_counters = []
        self.full_names_counters = []

        try:
            full_path_for_header_file = os.path.join(baseName, header_file)
        except:
            raise TypeError("You did not specify a directory path.")
        self.REIXSconfig = REIXSconfig

        # Open the header file and read line by line
        with open(full_path_for_header_file) as f:
            for line in f:
                if line.startswith('#S '):  # Triggers a new scan
                    # Append new scan if list is currently populated in last element (scan)
                    if self.datasets[-1] != []:
                        # we are in a new block
                        self.datasets.append([])

                    # Ensures we can only load scans with the same scan number only once
                    if line.strip().split()[1] in self.scanNumbers:
                        raise ValueError(
                            "Recorded Scans with the same scan number")

                    # Append scan information
                    else:
                        self.scanNumbers.append(line.strip().split()[1])
                        self.scanType.append(line.strip().split()[2])
                        self.scanMotor.append(line.strip().split()[3])

                # Header
                elif line.startswith("#L"):
                    self.datasets[-1].append(line.strip("#L ").strip("\n"))

                # Use this to generate mnemonic/long name dict
                elif line.startswith('#J'):
                    plist0 = line.strip().split(" ", 1)
                    plist = plist0[1].split("  ")
                    self.full_names_counters += plist

                elif line.startswith('#O'):
                    plist0 = line.strip().split(" ", 1)
                    plist = plist0[1].split("  ")
                    self.full_names_motors += plist

                # Use this to generate mnemonic/long name dict
                elif line.startswith('#j'):
                    plist0 = line.strip().split(" ", 1)
                    plist = plist0[1].split(" ")
                    self.mnemonics_counters += plist

                elif line.startswith('#o'):
                    plist0 = line.strip().split(" ", 1)
                    plist = plist0[1].split(" ")
                    self.mnemonics_motors += plist

                # Ignore empty spaces or commented lines
                elif line.startswith('#') or line.startswith('\n'):
                    pass

                # This is the actual data
                else:
                    self.datasets[-1].append(line.strip("\n"))

        # Now load the sdd data from file
        self.sdd_datasets = [[]]
        self.sdd_scanNumbers = []
        try:
            with open(f"{os.path.join(baseName,header_file)}"+"_sdd") as f_sdd:
                next(f_sdd)
                # Read file line by line
                for line_sdd in f_sdd:
                    if line_sdd.startswith('#S '):  # Triggers new scan
                        if self.sdd_datasets[-1] != []:
                            # we are in a new block
                            self.sdd_datasets.append([])
                        self.sdd_scanNumbers.append(
                            line_sdd.strip().split()[1])

                    elif line_sdd.startswith('#') or line_sdd.startswith('\n'):
                        pass

                    else:
                        self.sdd_datasets[-1].append(line_sdd)

            # Since we loaded data, set SDD data flag to true
            self.SDD = True
        except:
            self.SDD = False
            warnings.warn("No SDD file found.")

        # Load MCP data
        self.mcp_datasets = [[]]
        self.mcp_scanNumbers = []
        try:
            # Try opening MCP file for RIXS endstation
            with open(f"{os.path.join(baseName,header_file)}"+"_mcpMCA") as f_mcp:
                next(f_mcp)
                # Read file line by line
                for line_mcp in f_mcp:
                    if line_mcp.startswith('#S '):
                        if self.mcp_datasets[-1] != []:
                            # we are in a new block
                            self.mcp_datasets.append([])
                        self.mcp_scanNumbers.append(
                            line_mcp.strip().split()[1])

                    elif line_mcp.startswith('#') or line_mcp.startswith('\n'):
                        pass

                    else:
                        self.mcp_datasets[-1].append(line_mcp)

            # Set MCP data for RIXS endstation
            self.MCPRIXS = True
            self.MCPRSXS = False
        except:
            self.MCPRIXS = False
            try:
                # Try loading MCP data file from RSXS endstation
                self.mcpRSXS_data = [[]]
                self.mcpRSXS_scale_headers = [[]]
                self.mcpRSXS_scanNumbers = []

                with open(f"{os.path.join(baseName,header_file)}"+"_mcp") as f_mcp:
                    next(f_mcp)
                    # Read file line by line
                    for line in f_mcp:
                        if line.startswith('#S '):  # Triggers new scan
                            if self.mcpRSXS_data[-1] != []:
                                self.mcpRSXS_data.append([])
                            # Now need to also keep track of column headers
                            if self.mcpRSXS_scale_headers[-1] != []:
                                self.mcpRSXS_scale_headers.append([])

                            self.mcpRSXS_scanNumbers.append(
                                line.strip().split()[1])

                        # Append the header
                        elif line.startswith('#C ') or line.startswith('#@IMG'):
                            self.mcpRSXS_data[-1].append([])
                            if self.mcpRSXS_scale_headers[-1] == []:
                                for i in (line.replace("\t", "").strip("\n").split(' ')[1:3]):
                                    self.mcpRSXS_scale_headers[-1].append(i)

                        # Skip empty spaces
                        elif line.startswith('#') or line.startswith('\n'):
                            pass

                        # Append the actual data
                        else:
                            self.mcpRSXS_data[-1][-1].append(line.strip('\n'))
                self.MCPRSXS = True
            except:
                self.MCPRSXS = False
                warnings.warn("No MCP file found.")

        # Load XEOL data if it exists
        self.xeol_datasets = [[]]
        self.xeol_scanNumbers = []
        try:
            with open(f"{os.path.join(baseName,header_file)}"+"_xeol") as f_xeol:
                next(f_xeol)
                # Load file line by line
                for line_xeol in f_xeol:
                    if line_xeol.startswith('#S '):  # Triggers new scan
                        if self.xeol_datasets[-1] != []:
                            # we are in a new block
                            self.xeol_datasets.append([])
                        self.xeol_scanNumbers.append(
                            line_xeol.strip().split()[1])

                    # Ignore empty spaces
                    elif line_xeol.startswith('#') or line_xeol.startswith('\n'):
                        pass

                    # Append actual data
                    else:
                        self.xeol_datasets[-1].append(line_xeol)

            self.XEOL = True

        except:
            self.XEOL = False

        # Check for scan number mismatches across all ASCII files
        if self.scanNumbers != self.sdd_scanNumbers:
            print("Scan Number Mismatch in Files (SDD)")

        if self.MCPRIXS == True and self.scanNumbers != self.mcp_scanNumbers:
            print("Scan Number Mismatch in Files (RIXS MCP)")

        if self.MCPRSXS == True and self.scanNumbers != self.mcpRSXS_scanNumbers:
            print("Scan Number Mismatch in Files (RSXS MCP)")

        if self.XEOL == True and self.scanNumbers != self.xeol_scanNumbers:
            print("Scan Number Mismatch in Files (XEOL)")

        if int(self.scanNumbers[-1]) != int(len(self.scanNumbers)):
            print("Scans not labelled consecutively!")
            print("Last Scan Number", self.scanNumbers[-1])
            print("Scans in File", len(self.scanNumbers))

        if len(self.mnemonics_motors) != len(self.full_names_motors):
            raise UserWarning(
                "Mismatch with nmemonic to full name dictionary (motors).")

        if len(self.mnemonics_counters) != len(self.full_names_counters):
            raise UserWarning(
                "Mismatch with nmemonic to full name dictionary (counters).")

        self.mnemonic2name_motors = dict(
            zip(self.full_names_motors, self.mnemonics_motors))
        self.mnemonic2name_counters = dict(
            zip(self.full_names_counters, self.mnemonics_counters))

        # Create dictionary for scan numbers
        # This is in case a detector was disabled for a specific scan,
        # so we know which list index to trigger
        self.scanIndexDictHeader = dict()
        for i, scanIndex in enumerate(self.scanNumbers):
            self.scanIndexDictHeader[int(scanIndex)] = i

        if self.SDD == True:
            self.scanIndexDictSDD = dict()
            for i, scanIndex in enumerate(self.sdd_scanNumbers):
                self.scanIndexDictSDD[int(scanIndex)] = i

        if self.MCPRIXS == True:
            self.scanIndexDictRIXSMCP = dict()
            for i, scanIndex in enumerate(self.mcp_scanNumbers):
                self.scanIndexDictRIXSMCP[int(scanIndex)] = i

        if self.MCPRSXS == True:
            self.scanIndexDictRSXSMCP = dict()
            for i, scanIndex in enumerate(self.mcpRSXS_scanNumbers):
                self.scanIndexDictRSXSMCP[int(scanIndex)] = i

        if self.XEOL == True:
            self.scanIndexDictXEOL = dict()
            for i, scanIndex in enumerate(self.xeol_scanNumbers):
                self.scanIndexDictXEOL[int(scanIndex)] = i

    def Scan(self, scan):
        """Returns all data associated with a specific Scan.

        Parameters
        ----------
        scan : int
            The scan to be loaded

        Returns
        -------
        REIXS : An scan object that can be used for analysis
        """

        return self.myScan(self, scan)

    class myScan(object):
        """Returns all data associated with a specific Scan."""

        def __init__(my, self, scan):
            """
            Constructor for a specific Scan.

            Parameters
            ----------
            scan : int
                The scan to be loaded

            Returns
            -------
            myScan : An scan object that can be used for analysis
            """

            # Get index of requested scan
            # Do this for all data files
            try:
                my.scan = scan
                my.scanl = self.scanIndexDictHeader[scan]
            except:
                raise ValueError("Scan not defined in header file.")

            try:
                my.scansdd = self.scanIndexDictSDD[scan]
                my.SDD = True
            except:
                my.SDD = False

            try:
                my.scanmcpRIXS = self.scanIndexDictRIXSMCP[scan]
                my.MCPRIXS = True
            except:
                my.MCPRIXS = False

            try:
                my.scanmcpRSXS = self.scanIndexDictRSXSMCP[scan]
                my.MCPRSXS = True
            except:
                my.MCPRSXS = False

            try:
                my.scanxeol = self.scanIndexDictXEOL[scan]
                my.XEOL = True
            except:
                my.XEOL = False

            # Load header file / SCA data into pandas data frame
            # Need to load data and header separately, then combine
            my.sca_data = pd.DataFrame(
                self.datasets[my.scanl][1:]).iloc[:, 0].str.split(" ", expand=True)
            header = pd.DataFrame(self.datasets[my.scanl]).iloc[0].str.split(
                "  ", expand=True).transpose()
            header_list = [i for i in header[0]]
            epoch_index = header_list.index('Epoch')

            mnemonics_header = list()
            # Do the motors first
            for entry in header_list[0:epoch_index]:
                mnemonics_header.append(self.mnemonic2name_motors[entry])
            # Now append Epoch
            mnemonics_header.append("Epoch")
            # Do the counters last
            for entry in header_list[epoch_index+1:]:
                mnemonics_header.append(self.mnemonic2name_counters[entry])

            if len(header_list) != len(mnemonics_header):
                warnings.warn(
                    "Problem in converting long names to mnemonics in header file.")

            my.sca_data.columns = mnemonics_header
            # Ensures data type integrity
            my.sca_data = my.sca_data.apply(pd.to_numeric, errors='coerce')

            # Remove duplicate columns
            my.sca_data = my.sca_data.loc[:, ~
                                          my.sca_data.columns.duplicated()].copy()

            # Defines special streams
            # Added legacy support for previous spec variables
            try:
                my.mono_energy = np.array(
                    my.sca_data[self.REIXSconfig["ASCII_mono_energy"]])
            except:
                # We leave this in for legacy support
                try:
                    my.mono_energy = np.array(my.sca_data["beam"])
                except:
                    raise TypeError("Problem determining energy.")

            try:
                my.mesh_current = np.array(
                    my.sca_data[self.REIXSconfig["ASCII_mesh_current"]])
            except:
                # Leave this in for legacy support
                try:
                    my.mesh_current = np.array(my.sca_data["i0"])
                except:
                    raise TypeError("Problem determening mesh current")

            try:
                my.sample_current = np.array(
                    my.sca_data[self.REIXSconfig["ASCII_sample_current"]])
            except:
                # Also leave this in for legacy support
                try:
                    my.sample_current = np.array(my.sca_data["tey"])
                except:
                    raise TypeError("Problem determening sample current.")

            # Define total electron yield normalized by flux
            try:
                my.TEY = my.sample_current/my.mesh_current
            except:
                print("Error storing SCAs to pre-defined variables. Likely column headers have changed. You can access SCAs by column name.")

            # Load XAS / PFY data
            if my.SDD == True:
                try:
                    # Load the SDD MCA
                    my.sdd_data = np.loadtxt(
                        self.sdd_datasets[my.scansdd], skiprows=1024, dtype='float')
                    # Load the SDD energy scale
                    my.sdd_energy = np.loadtxt(
                        self.sdd_datasets[my.scansdd], max_rows=1024)
                except:
                    raise UserWarning("Could not load SDD data from file.")

            # Load the XES / MCP data for RIXS endstation
            if my.MCPRIXS == True:
                try:
                    # Load the MCP energy scale
                    my.mcp_energy = np.loadtxt(
                        self.mcp_datasets[my.scanmcpRIXS], max_rows=1024, dtype='float')
                    # Load the MCP MCA
                    my.mcp_data = np.loadtxt(
                        self.mcp_datasets[my.scanmcpRIXS], skiprows=1024, unpack=True)
                except:
                    raise UserWarning("Could not RIXS MCP data from file.")

            # Load the MCP data for RSXS endstation
            if my.MCPRSXS == True:
                try:
                    # Keep track of scales in MCP file
                    my.mcpRSXS_scales = dict()
                    my.mcpRSXS_scatters = dict()
                    my.mcpRSXS_axes = self.mcpRSXS_scale_headers[my.scanmcpRSXS]

                    j = 0
                    k = 0
                    # Scales and scatter images alternate in data file for one scan
                    # Append appropriately after separating
                    for i in range(len(self.mcpRSXS_data[my.scanmcpRSXS])):
                        if i % 2:  # is odd i
                            my.mcpRSXS_scatters[j] = np.loadtxt(
                                self.mcpRSXS_data[my.scanmcpRSXS][i])
                            j += 1

                        else:  # is even
                            my.mcpRSXS_scales[k] = np.transpose(
                                np.loadtxt(self.mcpRSXS_data[my.scanmcpRSXS][i]))
                            k += 1

                except:
                    raise UserWarning(
                        "Could not load RSXS MCP data from file.")

            # Load XEOL data
            if my.XEOL == True:
                try:
                    # Get the xeol energy on wavelength scale
                    my.xeol_energy = np.loadtxt(
                        self.xeol_datasets[my.scanxeol], max_rows=2048)
                    # Get the xeol background - this is not used (!)
                    my.xeol_background = np.loadtxt(
                        self.xeol_datasets[my.scanxeol], skiprows=2048, max_rows=2048)
                    # Load the actual xeol image
                    my.xeol_data = np.loadtxt(
                        self.xeol_datasets[my.scanxeol], skiprows=4096)
                except:
                    raise UserWarning("Could not load XEOL data from file.")
