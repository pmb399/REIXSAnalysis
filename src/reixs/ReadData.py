# Scientific Modules
import numpy as np
import pandas as pd

# Edge Dict
from .edges import EdgeDict
from .rsxs_readutil import img_to_sca,grid_stack,stack_to_mca

# Utilities
import os
import warnings

## Every scan will be an instance of a class ##

#########################################################################################

class REIXS(object):
    
    """Returns all scans in a scan file."""
    
    def __init__(self,baseName,header_file):
        """
        Constructor for one data file.
    
        Parameters:
            baseName: absolute path to experimental data.
            header_file: name of the header file.
            
        Returns:
            REIXS: all data associated with the specified data file ready to be processed.
    
        """
        
        ## Load all scans in specified file
        self.datasets = [[]]
        self.scanNumbers = []
        self.scanType = []
        self.scanMotor = []
        
        try:
            full_path_for_header_file = os.path.join(baseName,header_file)
        except:
            raise TypeError("You did not specify a directory path.")
            
        with open(full_path_for_header_file) as f:
            for line in f:
                if line.startswith('#S '):
                    if self.datasets[-1] != []:
                        # we are in a new block
                        self.datasets.append([])

                    if line.strip().split()[1] in self.scanNumbers:
                        raise ValueError("Recorded Scans with the same scan number")

                    else:
                        self.scanNumbers.append(line.strip().split()[1])
                        self.scanType.append(line.strip().split()[2])
                        self.scanMotor.append(line.strip().split()[3])

                elif line.startswith("#L"):
                    self.datasets[-1].append(line.strip("#L ").strip("\n"))

                elif line.startswith('#') or line.startswith('\n'):
                    pass

                else:
                    self.datasets[-1].append(line.strip("\n"))

        self.sdd_datasets = [[]]
        self.sdd_scanNumbers = []
        try:
            with open(f"{os.path.join(baseName,header_file)}"+"_sdd") as f_sdd:
                next(f_sdd)
                for line_sdd in f_sdd:
                    if line_sdd.startswith('#S '):
                        if self.sdd_datasets[-1] != []:
                            # we are in a new block
                            self.sdd_datasets.append([])
                        self.sdd_scanNumbers.append(line_sdd.strip().split()[1])

                    elif line_sdd.startswith('#') or line_sdd.startswith('\n'):
                        pass

                    else:
                        self.sdd_datasets[-1].append(line_sdd)

            self.SDD = True
        except:
            self.SDD = False
            warnings.warn("No SDD file found.")

        self.mcp_datasets = [[]]
        self.mcp_scanNumbers = []
        try:
            with open(f"{os.path.join(baseName,header_file)}"+"_mcpMCA") as f_mcp:
                next(f_mcp)
                for line_mcp in f_mcp:          
                    if line_mcp.startswith('#S '):
                        if self.mcp_datasets[-1] != []:
                            # we are in a new block
                            self.mcp_datasets.append([])
                        self.mcp_scanNumbers.append(line_mcp.strip().split()[1])

                    elif line_mcp.startswith('#') or line_mcp.startswith('\n'):
                        pass

                    else:
                        self.mcp_datasets[-1].append(line_mcp)

            self.MCPRIXS = True
            self.MCPRSXS = False
        except:
            self.MCPRIXS = False
            try:
                self.mcpRSXS_data = [[]]
                self.mcpRSXS_scale_headers = [[]]
                self.mcpRSXS_scanNumbers = []

                with open(f"{os.path.join(baseName,header_file)}"+"_mcp") as f_mcp:
                    next(f_mcp)
                    for line in f_mcp:
                        if line.startswith('#S '):
                            if self.mcpRSXS_data[-1]!=[]:
                                self.mcpRSXS_data.append([])
                            if self.mcpRSXS_scale_headers[-1] != []:
                                self.mcpRSXS_scale_headers.append([])
                            
                            self.mcpRSXS_scanNumbers.append(line.strip().split()[1])
                            
                        elif line.startswith('#C ') or line.startswith('#@IMG'):
                            self.mcpRSXS_data[-1].append([])
                            if self.mcpRSXS_scale_headers[-1] == []:
                                for i in (line.replace("\t","").strip("\n").split(' ')[1:3]):
                                    self.mcpRSXS_scale_headers[-1].append(i)
                                

                        elif line.startswith('#') or line.startswith('\n'):
                            pass

                        else:
                            self.mcpRSXS_data[-1][-1].append(line.strip('\n'))
                self.MCPRSXS = True
            except:
                self.MCPRSXS = False
                warnings.warn("No MCP file found.")

                    
        self.xeol_datasets = [[]]
        self.xeol_scanNumbers = []
        try:
            with open(f"{os.path.join(baseName,header_file)}"+"_xeol") as f_xeol:
                next(f_xeol)
                for line_xeol in f_xeol:          
                    if line_xeol.startswith('#S '):
                        if self.xeol_datasets[-1] != []:
                            # we are in a new block
                            self.xeol_datasets.append([])
                        self.xeol_scanNumbers.append(line_xeol.strip().split()[1])

                    elif line_xeol.startswith('#') or line_xeol.startswith('\n'):
                        pass

                    else:
                        self.xeol_datasets[-1].append(line_xeol)
                        
            self.XEOL = True
            #print("XEOL data found.")
            
            if self.scanNumbers != self.xeol_scanNumbers:
                print("Scan Number Mismatch in Files")

        except:
            self.XEOL = False


        if self.scanNumbers != self.sdd_scanNumbers:
            print("Scan Number Mismatch in Files (Header)")

        if self.MCPRIXS == True and self.scanNumbers != self.mcp_scanNumbers:
            print("Scan Number Mismatch in Files (RIXS MCP)")

        if self.MCPRSXS==True and self.scanNumbers != self.mcpRSXS_scanNumbers:
            print("Scan Number Mismatch in Files (RSXS MCP)")

        if self.XEOL == True and self.scanNumbers != self.xeol_scanNumbers:
            print("Scan Number Mismatch in Files (XEOL)")

        if int(self.scanNumbers[-1]) != int(len(self.scanNumbers)):
            print("Scans not labelled consecutively!")
            print("Last Scan Number",self.scanNumbers[-1])
            print("Scans in File", len(self.scanNumbers))
            
        ## Create dictionary for scan numbers
        self.scanIndexDictHeader = dict()
        for i,scanIndex in enumerate(self.scanNumbers):
            self.scanIndexDictHeader[int(scanIndex)] = i

        if self.SDD == True:
            self.scanIndexDictSDD = dict()
            for i,scanIndex in enumerate(self.sdd_scanNumbers):
                self.scanIndexDictSDD[int(scanIndex)] = i

        if self.MCPRIXS == True:
            self.scanIndexDictRIXSMCP = dict()
            for i,scanIndex in enumerate(self.mcp_scanNumbers):
                self.scanIndexDictRIXSMCP[int(scanIndex)] = i

        if self.MCPRSXS == True:
            self.scanIndexDictRSXSMCP = dict()
            for i,scanIndex in enumerate(self.mcpRSXS_scanNumbers):
                self.scanIndexDictRSXSMCP[int(scanIndex)] = i

        if self.XEOL == True:
            self.scanIndexDictXEOL = dict()
            for i,scanIndex in enumerate(self.xeol_scanNumbers):
                self.scanIndexDictXEOL[int(scanIndex)] = i

    def Scan(self,scan):
        """Returns all data associated with a specific Scan.
        
        Parameters:
            scan (REIXS) : The scan to be loaded
            
        Returns:
            REIXS : An scan object that can be used for analysis
        """        
        
        return self.myScan(self,scan)               
        
    class myScan(object):
        """Returns all data associated with a specific Scan.""" 

        def __init__(my,self,scan):
            """
            Constructor for a specific Scan.
        
            Parameters:
                scan (REIXS) : The scan to be loaded
            
            Returns:
                myScan : An scan object that can be used for analysis
            """ 
            
            try:
                my.scan  = scan
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

            ## Load header file / SCA data
            my.sca_data = pd.DataFrame(self.datasets[my.scanl][1:]).iloc[:,0].str.split(" ", expand=True)
            header  = pd.DataFrame(self.datasets[my.scanl]).iloc[0].str.split("  ", expand=True).transpose()
            my.sca_data.columns = header[0]
            my.sca_data = my.sca_data.apply(pd.to_numeric, errors='coerce')
                
            try:
                my.mono_energy = np.array(my.sca_data["Mono_Engy"])
            except:
                try:
                    my.mono_energy = np.array(my.sca_data["Beam Engy"])
                except:
                    try:
                        my.mono_energy = np.array(my.sca_data["Mono Ener"])
                    except:
                        try:
                            my.mono_energy = np.array(my.sca_data["BeamEngy"])
                        except:
                            raise TypeError("Problem determining energy.")

                
            try:
                my.mesh_current = np.array(my.sca_data["Mesh"])
            except:
                try:
                    my.mesh_current = np.array(my.sca_data["Mesh Curr"])
                except:
                    try:
                        my.mesh_current = np.array(my.sca_data["I0_BD3"])
                    except:
                        raise TypeError("Problem determening mesh current")

            try:
                my.sample_current = np.array(my.sca_data["Sample"])
            except:
                try:
                    my.sample_current = np.array(my.sca_data["Samp Curr"])
                except:
                    try:
                        my.sample_current = np.array(my.sca_data["TEY"])
                    except:
                        raise TypeError("Problem determening sample current.")

            try:
                my.TEY = my.sample_current/my.mesh_current
            except:
                print("Error storing SCAs to pre-defined variables. Likely column headers have changed. You can access SCAs by column name.")


            ## Load XAS / PFY data
            if my.SDD == True:
                try:
                    my.sdd_data = np.loadtxt(self.sdd_datasets[my.scansdd],skiprows=1024,dtype='float')
                    my.sdd_energy = np.loadtxt(self.sdd_datasets[my.scansdd], max_rows=1024)
                except:
                    raise UserWarning("Could not load SDD data from file.")

            ## Load the XES / MCP data
            if my.MCPRIXS == True:
                try:
                    my.mcp_energy = np.loadtxt(self.mcp_datasets[my.scanmcpRIXS], max_rows=1024, dtype='float')
                    my.mcp_data = np.loadtxt(self.mcp_datasets[my.scanmcpRIXS],skiprows=1024,unpack=True)
                except:
                    raise UserWarning("Could not RIXS MCP data from file.")

            if my.MCPRSXS == True:
                try:
                    my.mcpRSXS_scales = dict()
                    my.mcpRSXS_scatters = dict()
                    my.mcpRSXS_axes = self.mcpRSXS_scale_headers[my.scanmcpRSXS]

                    j = 0
                    k = 0
                    for i in range(len(self.mcpRSXS_data[my.scanmcpRSXS])):
                        if i % 2: # is odd i
                            my.mcpRSXS_scatters[j] = np.loadtxt(self.mcpRSXS_data[my.scanmcpRSXS][i])
                            j += 1
                            
                        else: # is even
                            my.mcpRSXS_scales[k] = np.transpose(np.loadtxt(self.mcpRSXS_data[my.scanmcpRSXS][i]))
                            k += 1

                except:
                    raise UserWarning("Could not load RSXS MCP data from file.")
            
            ## Load XEOL data
            if my.XEOL==True:
                try:
                    my.xeol_energy = np.loadtxt(self.xeol_datasets[my.scanxeol],max_rows=2048)
                    my.xeol_background = np.loadtxt(self.xeol_datasets[my.scanxeol],skiprows=2048,max_rows=2048)
                    my.xeol_data = np.loadtxt(self.xeol_datasets[my.scanxeol],skiprows=4096)
                except:
                    raise UserWarning("Could not load XEOL data from file.")
                

        def MCP_norm(my):
            """Normalize the counts of the MCP by incident flux at every given datapoint.
               This is only applied to EEMs."""
            mcp_data_norm_init = []
            for count,i in enumerate(my.mesh_current):
                mcp_data_norm_i = np.true_divide(my.mcp_data[:,count],i)
                mcp_data_norm_init.append(mcp_data_norm_i)
            my.mcp_data_norm = (np.array(mcp_data_norm_init))

            return my.mcp_data_norm

        def XES(my):
            """Sum the MCP detector image over all recorded datapoints."""
            return np.sum(my.mcp_data,axis=1)
        
        def rXES(my,xes_incident_start,xes_incident_end):
            """Calculate resonant emission at selected energy (MCP)."""
            idx_start = (np.abs(xes_incident_start - my.mono_energy)).argmin()
            idx_end = (np.abs(xes_incident_end - my.mono_energy)).argmin()

            if not(hasattr(my, 'mcp_data_norm')):
                my.MCP_norm()

            return my.mcp_data_norm[idx_start:idx_end,:].sum(axis=0)
        
        def SDD_norm(my):
            sdd_data_norm_init = list()
            for i,mesh_current in enumerate(my.mesh_current):
                sdd_data_norm_init.append(np.true_divide(my.sdd_data[i,:],mesh_current))
            my.sdd_data_norm = np.array(sdd_data_norm_init)
            
            return my.sdd_data_norm
        
        def XRF(my):
            """Sum the SDD detector image over all recorded datapoints"""
            return np.sum(my.sdd_data,axis=0)
        
        def rXRF(my,xes_incident_start,xes_incident_end):
            """Calculate resonant emission at selected energy (SDD)."""
            idx_start = (np.abs(xes_incident_start - my.mono_energy)).argmin()
            idx_end = (np.abs(xes_incident_end - my.mono_energy)).argmin()

            if not(hasattr(my, 'sdd_data_norm')):
                my.SDD_norm()

            return my.sdd_data_norm[idx_start:idx_end,:].sum(axis=0)
        
#         def XEOL_remove_bg(my):
#             my.xeol_data_bgrm = np.copy(my.xeol_data)
            
#             for i,arr in enumerate(my.xeol_data):
#                 my.xeol_data_bgrm[i,:] = np.subtract(arr,my.xeol_background)
                
#             return my.xeol_data_bgrm
        
#         def XEOL_norm_counts(my,bgrm):
#             if bgrm == True:
#                 if not hasattr(my,"xeol_data_bgrm"):
#                     my.XEOL_remove_bg()
#                 my.xeol_data_norm = np.copy(my.xeol_data_bgrm)
#             else:
#                 my.xeol_data_norm = np.copy(my.xeol_data)
            
#             for i,counts in enumerate(my.mesh_current):
#                 my.xeol_data_norm[i,:] = np.true_divide(np.max(my.mesh_current)*my.xeol_data[i,:],counts)
                
#             return my.xeol_data_norm
        
#         def XEOL_spec(my,bgrm,fluxnorm):
#             if bgrm == True:
#                 if fluxnorm == True:
#                     xeol_data = my.XEOL_norm_counts(bgrm=True)
#                 else:
#                     xeol_data = my.XEOL_remove_bg()
#             else:
#                 if fluxnorm == True:
#                     xeol_data = my.XEOL_norm_counts(bgrm=False)
#                 else:
#                     xeol_data = my.xeol_data
            
#             return np.sum(xeol_data,axis=0)

        
        def ROI(my,PFY_edge=None,iPFY_edge=None):
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
                
                return SDDLowerBound,SDDUpperBound
                
            elif iPFY_edge != None:
                iSDDLowerBound = EdgeDict[iPFY_edge][0]
                iSDDUpperBound = EdgeDict[iPFY_edge][1]
                
                return iSDDLowerBound,iSDDUpperBound
            
            else:
                raise Exception("Proper ROI bounds not defined.")
                

        def PFY(my,PFY_edge=None,SDDLowerBound=None,SDDUpperBound=None):
            """ Calculate PFY based on a defined ROI.
            
            Prerequisites:
                Option 1: Define ROI via ROI method
                Option 2: Define ROI explictily by calling
                            SDDLowerBound : lower energy limit for summation
                            SDDUpperBound : upper energy limit for summation
            
            """
            if PFY_edge != None:
                SDDLowerBound,SDDUpperBound = my.ROI(PFY_edge)
                
            elif SDDLowerBound != None and SDDUpperBound != None:
                pass
            
            else:
                raise Exception("Proper PFY boundaries not defined. Define boundaries or use TFY.")
                
            idx_low_sdd = (np.abs(SDDLowerBound - my.sdd_energy)).argmin()
            idx_high_sdd = (np.abs(SDDUpperBound - my.sdd_energy)).argmin()

            ## Transform this to PFY spectrum (rows = incident energy points, columns = detector 'space')
            sdd_detector_sum = my.sdd_data[:,idx_low_sdd:idx_high_sdd].sum(axis=1)
            PFY_spec = (sdd_detector_sum/my.mesh_current)

            return PFY_spec

        def iPFY(my,iPFY_edge=None,iSDDLowerBound=None,iSDDUpperBound=None):
            """ Calculate iPFY based on a defined ROI.
            
                Prerequisites:
                    Option 1: Define ROI via ROI method (and specify **kwarg iPFY_edge)
                    Option 2: Define ROI explictily by calling
                            iSDDLowerBound : lower energy limit for summation
                            iSDDUpperBound : upper energy limit for summation
            
            """
            if iPFY_edge != None:
                iSDDLowerBound,iSDDUpperBound = my.ROI(iPFY_edge)
                
            elif iSDDLowerBound != None and iSDDUpperBound != None:
                pass
            
            else:
                raise Exception("Proper iPFY boundaries not defined.")

            iidx_low_sdd =  (np.abs(iSDDLowerBound - my.sdd_energy)).argmin()
            iidx_high_sdd = (np.abs(iSDDUpperBound  - my.sdd_energy)).argmin()
            sdd_detector_sum_i = my.sdd_data[:,iidx_low_sdd:iidx_high_sdd].sum(axis=1)
            iPFY_spec = 1/(sdd_detector_sum_i/my.mesh_current)

            return iPFY_spec
            
        def TFY(my):
            """Calculate TFY by summing over entire SDD image without ROI."""
            ## Transform this to TFY spectrum (rows = incident energy points, columns = detector 'space')
            sdd_detector_sum = my.sdd_data[:,:].sum(axis=1)
            TFY_spec = (sdd_detector_sum/my.mesh_current)

            return TFY_spec
                

        def specPFY(my,mcp_lowE,mcp_highE):
            """Calculate spectrometer PFY based on ROI set."""
            mcp_highE_idx = (np.abs(mcp_highE-my.mcp_energy)).argmin()
            mcp_lowE_idx = (np.abs(mcp_lowE-my.mcp_energy)).argmin()

            if not(hasattr(my, 'mcp_data_norm')):
                my.MCP_norm()

            return my.mcp_data_norm[:,mcp_lowE_idx:mcp_highE_idx].sum(axis=1)

        def RSXS_MCPnorm(my):
            my.RSXSMCP_norm = dict()
            for k,v in my.mcpRSXS_scatters.items():
                my.RSXSMCP_norm[k] = np.true_divide(v,my.sample_current[k])

            return my.RSXSMCP_norm
                        
        def RSXS_1dROI(my,img,x_low=None,x_high=None,y_low=None,y_high=None,axis=0):

            if not(hasattr(my, 'RSXSMCP_norm')):
                my.RSXS_MCPnorm()

            return img_to_sca(my.mcpRSXS_scales,my.RSXSMCP_norm,img,x_low,x_high,y_low,y_high,axis)

        def RSXS_2dROI(my,x_low=None,x_high=None,y_low=None,y_high=None,axis=0):
            if not(hasattr(my, 'RSXSMCP_norm')):
                my.RSXS_MCPnorm()

            return stack_to_mca(my.mcpRSXS_scales,my.RSXSMCP_norm,x_low,x_high,y_low,y_high,axis)

        def RSXS_Images(my):
            if not(hasattr(my, 'RSXSMCP_norm')):
                my.RSXS_MCPnorm()

            return grid_stack(my.mcpRSXS_scales,my.RSXSMCP_norm)