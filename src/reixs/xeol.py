import numpy as np

#########################################################################################
def XEOL_background_removal(data,arg,REIXSobj,background_scan=None):
    
    if background_scan!=None:
        # 1 Sum the background
        background_data = REIXSobj.Scan(background_scan).xeol_data
        background_sum = background_data.sum(axis=0)

        # 2 Normalize to average background frame
        background_spec = np.true_divide(background_sum,len(background_data))

        # 3 Subtract the background from the actual spectrum
        xeol_data = data[arg].xeol_data
        # To do so, need to broadcast background into shape
        xeol_subtracted = np.subtract(xeol_data,background_spec)

        return xeol_subtracted

    else:
        return data[arg].xeol_data

def xeol_spec(data,arg,REIXSobj,background_scan=None):
    xeol_subtracted = XEOL_background_removal(data,arg,REIXSobj,background_scan)
    return np.sum(xeol_subtracted,axis=0)

def xeol_idx_lambda(data,arg,lambda_low,lambda_high):
    idx_low  = (np.abs(lambda_low  - data[arg].xeol_energy)).argmin()
    idx_high = (np.abs(lambda_high - data[arg].xeol_energy)).argmin()
    return idx_low,idx_high

        
def xeol_idx_energy(data,arg,energy_low,energy_high):
    idx_low  = (np.abs(energy_low  - data[arg].mono_energy)).argmin()
    idx_high = (np.abs(energy_high - data[arg].mono_energy)).argmin()
    return idx_low,idx_high


def rxeol_spec(data,arg,REIXSobj,energy_low,energy_high,background_scan=None):
    xeol_subtracted = XEOL_background_removal(data,arg,REIXSobj,background_scan)
    ixd_low,idx_high = xeol_idx_energy(data,arg,energy_low,energy_high)
    return np.sum(np.true_divide(xeol_subtracted[ixd_low:idx_high,:],data[arg].mesh_current[ixd_low:idx_high,None]),axis=0)

def toy_spec(data,arg,REIXSobj,background_scan=None):
    xeol_subtracted = XEOL_background_removal(data,arg,REIXSobj,background_scan)
    return np.sum(np.true_divide(xeol_subtracted,data[arg].mesh_current[:,None]),axis=1)

def poy_spec(data,arg,REIXSobj,lambda_low,lambda_high,background_scan=None):
    xeol_subtracted = XEOL_background_removal(data,arg,REIXSobj,background_scan)
    ixd_low,idx_high = xeol_idx_lambda(data,arg,lambda_low,lambda_high)
    return np.sum(np.true_divide(xeol_subtracted[:,ixd_low:idx_high],data[arg].mesh_current[:,None]),axis=1)


def xeol_map(data,arg,REIXSobj,background_scan=None):
    xeol_subtracted = XEOL_background_removal(data,arg,REIXSobj,background_scan)
    return np.true_divide(xeol_subtracted,data[arg].mesh_current[:,None])
