import os
import ast

def get_REIXSconfig():
    """Allows to read in an external spec config, otherwise default to internal value below
    
        Instructions
        ------------
        0) Import os
        1) Set up a REIXSconfig dictionary
            REIXSconfig = dict() 
        2) Append all python variables as keys into dict with value as the spec variabe
            REIXSconfig["python_variable"] = "spec variable"
        3) Export the REIXSconfig dict to environment
            os.environ["REIXSconfig"] = str(REIXSconfig)
    """
    # Read in external configuration if available
    if "REIXSconfig" in os.environ:
        # Expects environmental variable as string - need to convert to dict
        return ast.literal_eval(os.environ["REIXSconfig"])
    else: # default to internal values
        REIXSconfig = dict()
        REIXSconfig["HDF5_mono_energy"] = "Data/beam"
        REIXSconfig["HDF5_mesh_current"] = "Data/i0"
        REIXSconfig["HDF5_sample_current"] = "Data/tey"
        REIXSconfig["HDF5_sdd_data"] = 'Data/sddMCA'
        REIXSconfig["HDF5_sdd_energy"] = 'Data/sddMCA_scale'
        REIXSconfig["HDF5_xeol_data"] = 'Endstation/Detectors/XEOL/xeolMCA'
        REIXSconfig["HDF5_xeol_energy"] = 'Endstation/Detectors/XEOL/xeolMCA_scale'
        REIXSconfig["HDF5_xeol_background"] = 'Endstation/Detectors/XEOL/xeolMCA_back'
        REIXSconfig["HDF5_mcp_data"] = 'Data/mcpMCA'
        REIXSconfig["HDF5_mcp_energy"] = 'Data/mcpMCA_scale'
        REIXSconfig["HDF5_sca_data"] = 'Data'

        REIXSconfig["ASCII_mono_energy"] = 'beam'
        REIXSconfig["ASCII_mesh_current"] = 'i0'
        REIXSconfig["ASCII_sample_current"] = 'tey'
        return REIXSconfig