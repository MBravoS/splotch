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
	
	from .defaults import Params
	from .base_func import val_checker
	from matplotlib.pyplot import style
	
	(path)
	
	### SPLOTCH specific
	# Get list of available parameters within the splotch.Params class
	validParams=[key for key in Params.__dict__.keys() if key[:1] != '_']
	
	# find first Splotch line
	begin=None
	with open(path, 'r') as file:
		for num, line in enumerate(file.readlines()):
			if ("##### splotch configuration" in line.lower()):
				begin=num
				
	spltFile=open(path, 'r')
	allLines=spltFile.readlines()
	spltLines=allLines[begin:]
	pltLines=allLines[:begin]
	spltFile.close()
	pltDict ={}
	
	for line in spltLines:
		if (":" in line): # This line calls a parameter
			par=line.split(":")[0].rstrip()
			if par[0] == '#':
				continue
			else:
				val=line.split(":")[-1].lstrip().replace("\n","")
				
				# Adjust parameter value for booleans and numbers (i.e. floats/integer)
				if ',' in val:
					val=[val_checker(v) for v in val.split(',')]
				else:
					val=val_checker(val)
				
				if (par in validParams): # only edit parameters that exist within Params Class
					setattr(Params,par,val)
	
	for line in pltLines:
		if (":" in line): # This line calls a parameter
			par=line.split(":")[0].rstrip().lstrip()
			if par[0] == '#':
				continue
			else:
				val=line.split('#')[0].split(":")[-1].lstrip().replace("\n","").rstrip()
				
				# Adjust parameter value for booleans and numbers (i.e. floats/integer)
				if ',' in val and "'" not in val:
					val=[val_checker(v.rstrip()) for v in val.split(',')]
				else:
					val=val_checker(val.rstrip())
				
				pltDict[par]=val
	style.use(pltDict)
	
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
	
	basedir=os.path.dirname(splotch.__file__)
	splotch.use_style("{0}/styles/default.style".format(basedir))
	
	return(None)

