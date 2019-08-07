def adjust_text(which=['x','y'],ax=None,text_kw={},**kwargs):

	"""Function which allows the user to adjust texts of one or many of either x/y-axis labels,
	title, legend, colorbars, etc.
	
	Parameters
	----------
	which : str, array-like or list
		Which `Text` instance(s) to apply the text properties to.

		Valid arguments can be one or many of the following:
			* 'x'|'xlabel'	 : x-axis label 
			* 'y'|'ylabel'	 : y-axis label 
			* 't'|'title'	 : Title
			* 's'|'suptitle' : Sup. title
			* 'l'|'legend'	 : Legend text
			* 'c'|'colorbar' : Color bar
			* 'T'|'text' 	 : Text objects
			* 'a'|'all'		 : All instances of all the above
	ax : pyplot.Axes or list, optional
		Use the given axes to adjust text, defaults to the current axis.
		If a list of Axis instances given, the text properties is applied to each.
	text_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib `Text` instance.
		It is recommended that text keyword arguments be given as **kwargs.
	**kwargs : Text instance properties
		kwargs are used to specify properties of `Text` instances
		A list of valid `Text` kwargs can be found here:
		[https://matplotlib.org/3.1.0/api/text_api.html#matplotlib.text.Text](https://matplotlib.org/3.1.0/api/text_api.html#matplotlib.text.Text "Matplotlib.text.Text")
	
	Returns
	-------
	None
	"""

	from matplotlib.text import Text
	from matplotlib.pyplot import gca
	from numpy import shape, min, max, argmin, argmax
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	try:
		_ = (it for it in ax)
	except TypeError:
		if (ax == None):
			ax = gca()

		ax = [ax]

	# Validate `which` value(s)
	whichRef = ['x','y','t','s','l','c','T','a',
				'xlabel','ylabel','title','suptitle','legend','colorbar','text','all']

	for w in which:
		try:
			wInd = whichRef.index(w)
			wComp = (wInd + len(whichRef)//2) % len(whichRef) # get the complimenting short/long version
			if (whichRef[wComp] in which):
				raise TypeError("adjust_text() received equivalent values for 'which': '{0}' and '{1}'.".format(whichRef[wInd],whichRef[wComp]))

		except (ValueError):
			if (type(which) != Text):
				raise TypeError("adjust_text() received invalid value for 'which' ('{0}'). Must be one of: {1}".format(w,', '.join(whichRef)))
		

	L = len(which)
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	textpar = text_kw.copy()
	textpar.update(kwargs)

	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	textpar = dict_splicer(textpar,L,[1]*L)

	for a in ax:
		for ii, lab in enumerate(which):
			if (lab in ['x','xlabel']):
				texts = [a.xaxis.label]
			elif (lab in ['y','ylabel']):
				texts = [a.yaxis.label]
			elif (lab in ['t','title']):
				texts = [a.title]
			elif (lab in ['s','suptitle']): # Not implemented
				texts = [a.title]
			elif (lab in ['l','legend']):
				texts = a.legend().get_texts() # Get list of text objects in legend
			elif (lab in ['c','colorbar']):
				# Get the axis with the largest ratio between width or height
				caxInd = argmax([max([c.get_position().width/c.get_position().height,c.get_position().height/c.get_position().width]) for c in a.figure.get_axes()])
				if (a.figure.get_axes()[caxInd].get_position().height > a.figure.get_axes()[caxInd].get_position().width):
					texts = [a.figure.get_axes()[caxInd].yaxis.label]
				else:
					texts = [a.figure.get_axes()[caxInd].xaxis.label]
			elif (lab in ['T','text']):
				texts = [child for child in a.get_children()[:-4] if type(child) == Text] # -4 to avoid grabbing Title, Subtitle. etc. Text instances
			elif (type(lab) == Text):
				texts = [lab]
			elif (lab in ['a','all']):
				texts = [a.xaxis.label,a.yaxis.label,a.title,a.legend().get_texts()]

				caxInd = argmax([max([c.get_position().width/c.get_position().height,c.get_position().height/c.get_position().width]) for c in a.figure.get_axes()])
				if (a.figure.get_axes()[caxInd].get_position().height > a.figure.get_axes()[caxInd].get_position().width):
					texts.append(a.figure.get_axes()[caxInd].yaxis.label)
				else:
					texts.append(a.figure.get_axes()[caxInd].xaxis.label)

				texts = texts + [child for child in a.get_children()[:-4] if type(child) == Text]

			for t in texts:
				t.set(**textpar[ii])
	
	return(None)


def subplots(naxes=None,nrows=None,ncols=None,va='top',ha='left',wspace=None,hspace=None,figsize=None,axes_kw={},**kwargs):
	
	from .base_func import dict_splicer
	
	from matplotlib import rcParams
	from matplotlib.gridspec import GridSpec
	from matplotlib.pyplot import figure, subplot
	
	from numpy import ceil
	
	gridRef = [[0,0], [1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [2,5], [3,4], [3,4], 
			   [4,4], [3,5], [3,5], [4,4], [3,6], [3,6], [4,5], [4,5], [3,7], [5,5], [5,5], [4,6], [5,5]]

	if (naxes == None): # No number of axes specified
		if (nrows==None): nrows=1
		if (ncols==None): ncols=1
		
		naxes = nrows*ncols
	else:
		if (ncols == None and nrows == None):
			if (naxes <= 25):
				nrows, ncols = gridRef[naxes] # Get best combination of rows x cols for number of axes
			else:
				raise ValueError(f"The naxes parameter is currently not implemented for naxes > 25.")
		elif (ncols != None and nrows != None):
			if (ncols*nrows != naxes):
				raise ValueError(f"Invalid number of axes ({naxes}) given for number of rows ({nrows}) and columns ({ncols}).") 
		else:
			if (nrows != None):
				ncols = int(ceil(naxes/nrows))
			else:
				nrows = int(ceil(naxes/ncols))
				
	
	# How many axes away from filling the gridspec evenly and completely
	delta = (nrows*ncols) - naxes
	
	if (figsize==None): # Auto scale default figure size to num cols/rows.
		figsize = (rcParams["figure.figsize"][0]*ncols, rcParams["figure.figsize"][1]*nrows)
	
	fig = figure(figsize=figsize,**kwargs)
	gs = GridSpec(ncols=ncols*2, nrows=nrows*2, hspace=hspace, wspace=wspace)
	
	axes_kw = dict_splicer(axes_kw,naxes,[1]*naxes)
	
	axes = []
	for ii in range(naxes):
		if (False):
			row = 2*(ncols-(ii%ncols)-1) if ha=='right' else 2*(ii%ncols)
			col = 2*(nrows-(ii//ncols)-1) if va=='bottom' else 2*(ii//ncols)
		else:
			row = 2*(nrows-(ii//ncols)-1) if va=='bottom' else 2*(ii//ncols)
			col = 2*(ncols-(ii%ncols)-1) if ha=='right' else 2*(ii%ncols)
			
		if (row == (0 if va=='bottom' else (nrows-1)*2) and ha=='centre'):
			axes.append(subplot(gs[row:row+2,col+delta:col+delta+2],**axes_kw[ii]))
		elif (col == (0 if ha=='right' else (ncols-1)*2) and va=='centre'):
			axes.append(subplot(gs[row+delta:row+delta+2,col:col+2],**axes_kw[ii]))
		else:
			axes.append(subplot(gs[row:row+2,col:col+2],**axes_kw[ii]))

	return(fig, axes)