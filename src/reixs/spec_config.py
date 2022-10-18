import os
import ast

def get_REIXSconfig():
    if "REIXSconfig" in os.environ:
        # Expects environmental variable as string - need to convert to dict
        return ast.literal_eval(os.environ["REIXSconfig"])
    else:
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
        REIXSconfig["HDF5_sca_data"] = 'Endstation/Counters'

        REIXSconfig["ASCII_mono_energy"] = 'Mono_Engy'
        REIXSconfig["ASCII_mesh_current"] = 'Mesh'
        REIXSconfig["ASCII_sample_current"] = 'Sample'
        return REIXSconfig