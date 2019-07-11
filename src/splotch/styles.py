def use_style(path):
    """ Loads in splotch configuration file containing a set of splotch style 
        parameters appended on top of a standard matplotlib .mplstyle configuration
        file.
    
    Parameters
    ----------
    path : string
            The path to the configuration file.

    Returns
    -------
    None
    
    """
    
    from matplotlib.pyplot import style
    from splotch.defaults import Params
    
    style.use(path)
    
    ### SPLOTCH specific
    # Get list of available parameters within the splotch.Params class
    validParams = [key for key in Params.__dict__.keys() if key[:1] != '_']
    
    # find first Splotch line
    begin = None
    with open(path, 'r') as file:
        for num, line in enumerate(file.readlines()):
            if ("##### splotch configuration" in line.lower()):
                begin = num
                
    spltFile = open(path, 'r')
    spltLines = spltFile.readlines()[begin:]
    
    for line in spltLines:
        if (":" in line): # This line calls a parameter
            par = line.split(":")[0].rstrip().replace("#","")
            val = line.split(":")[-1].lstrip().replace("\n","")
            
            # Adjust parameter value for booleans and numbers (i.e. floats/integer)
            try: # Check if val is a number
                float(val)
                isnum = True
            except ValueError:
                isnum = False
                
            if (val == "True"): val = True
            elif (val == "False"): val = False
            elif (isnum):
                val = float(val) if '.' in val else int(val)
            
            if (par in validParams): # only edit parameters that exist within Params Class
                setattr(Params,par,val)
    
    return(None)


def reset_style():
    """ Resets the Splotch parameters to their defaults.
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    
    """
    import os
    import splotch
    
    basedir = os.path.dirname(splotch.__file__)
    splotch.use_style("{0}/styles/default.style".format(basedir))
    
    return(None)
    
    
