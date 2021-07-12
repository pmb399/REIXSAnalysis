import numpy as np
import re
import parser
from reixs.ReadData import REIXS
from reixs.xeol import *
from reixs.offset import apply_offset

def loadMCAscans(basedir,file,x_stream,y_stream,detector,*args,norm=True,xoffset=None,xcoffset=None,yoffset=None,ycoffset=None,background=None):   
    ## Note that the data dict only gets populated locally until the appropriate
    ## singular y stream is return -- Assignment of y_stream happens after evaluation
    ## of mathemtatical expression
    
    def get_y_data(detector,data,arg,background,REIXSobj):
        ## Sets the detector
        if detector == "SDD Energy":
            y_data = data[arg].sdd_energy
        elif detector == 'MCP Energy':
            y_data = data[arg].mcp_energy
        elif detector == "XEOL Energy":
            y_data = data[arg].xeol_energy
            
        else:
            raise TypeError("Stream undefined!")

        return y_data

    def get_z_data(detector,data,arg,background,REIXSobj):
        ## Sets the detector
        if detector == "SDD":
            z_data = np.transpose(data[arg].SDD_norm())
        elif detector == 'MCP':
            z_data = np.transpose(data[arg].MCP_norm())
        elif detector == 'XEOL':
            if background!=None:
                z_data = np.transpose(xeol_map(data,arg,REIXSobj,background_scan=background))
            else:
                z_data = np.transpose(xeol_map(data,arg,REIXSobj))
            
        else:
            raise TypeError("Detector undefined!")

        return z_data

        
    def get_x_data(x_stream,data,arg,background,REIXSobj):
        if x_stream == "Mono Energy":
            return data[arg].mono_energy

        else:
            return np.array(data[arg].sca_data[x_stream])

        
    def math_stream(formula,data,arg,get_data,background=None,REIXSobj=None):
        ## Split the user input string at all mathematical operations
        ## Allow "( ) * / + -" as math
        
        pattern = '[\(+\-*^/\)]'
        split_expr = re.split(pattern,formula)

        ## Place string literals in dict if cannot be converted to float
        ## drop all empty strings from re.split
        quantity_str_dict = dict()

        for i,string in enumerate(split_expr):
            if string != "":
                try:
                    float(string)
                except:
                    # Assign generic "val{i}" key to string literal in compliance with
                    # python supported syntax for variables
                    quantity_str_dict[f"val{i}"] = string

        ## Parser does not support special string literals (ROIs) due to inproper python variable naming syntax
        ## Replace them with generic key as per dictionary
        ## Create local variable and assign corresponding data -- needed as eval interprets all input as variable
        ## This is indeed local to this function only and cannot be accessed from loadSCAscans
        for k,v in quantity_str_dict.items():
            formula = formula.replace(v,k)
            locals()[k] = get_data(v,data,arg,background,REIXSobj=REIXSobj)

        # Return the calculated result
        code = parser.expr(formula).compile()
        return eval(code)

    data = dict()
    REIXSobj = REIXS(basedir,file)
    for arg in args:
        data[arg] = REIXSobj.Scan(arg)

        ## Assign the calculated result to the y_stream of the object in data dict
        data[arg].y_data = math_stream(y_stream,data,arg,get_y_data)
        data[arg].y_data = apply_offset(data[arg].y_data,yoffset,ycoffset)
        
        ## Assign the calculated result to the x_stream of the object in data dict
        data[arg].x_data = math_stream(x_stream,data,arg,get_x_data)
        data[arg].x_data = apply_offset(data[arg].x_data,xoffset,xcoffset)

        data[arg].detector = math_stream(detector,data,arg,get_z_data,background,REIXSobj)

        if norm==True:
            data[arg].detector = data[arg].detector/np.max(data[arg].detector)
                
    return data
