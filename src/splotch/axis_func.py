def adjust_text(which=['x','y'],ax=None,text_kw={},**kwargs):

	""" Adjusts text instances.

	Function which allows the user to adjust texts of one or many of either x/y-axis labels,
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

	try: # check if iterable
		_ = (i for i in which)
		if (type(which) == str):
			which = [which]

	except (TypeError):
		which = [which]

	for w in which:
		try:
			wInd = whichRef.index(w)
			wComp = (wInd + len(whichRef)//2) % len(whichRef) # get the complimenting short/long version
			if (whichRef[wComp] in which):
				raise TypeError("adjust_text() received equivalent values for 'which': '{0}' and '{1}'.".format(whichRef[wInd],whichRef[wComp]))

		except (ValueError):
			if (type(w) != Text):
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


def subplots(naxes=None,nrows=None,ncols=None,va='top',ha='left',wspace=None,hspace=None,sharex='none',sharey='none',squeeze=True,figsize=None,axes_kw={},**kwargs):
	""" Adds a set of subplots to figure

	This is a more-generalised wrapper around matplotlib.pyplot.subplot function to allow for irregularly divided grids.
	
	Parameters
	----------
	naxes : int, optional, default: 1
		The number of axes objects to create.
		The resulting grid formed from specifying naxes is decided by ncols and nrows. The options are:

			- ncols and/or nrows not None:
				Makes sure that naxes can be correctly mapped into the specified grid.
				If one of nrows or ncols not given, the smallest possible grid will be made.
			- both ncols and nrows are None:
				Decides the best possible grid for this number of axes. Currently, this decision is hard-coded
				with plans for it to become an automatic decision later.
	
	nrows, ncols : int, optional
		Number of rows/columns of the subplot grid.
	
	va, ha : str, optional, default: 'top', 'left'
		The vertical alignment (va) and horizontal alignment (ha) sets the alignment of grids in the vertical and
		horizontal directions. 

		ha: 'left'		ha: 'centre'	ha: 'right'
		va: 'top'		va: 'top'		va: 'top'

		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯
		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯
		 ▯ ▯ ▯			  ▯ ▯ ▯			   ▯ ▯ ▯
		 ▯ ▯ ▯			  ▯ ▯ ▯			   ▯ ▯ ▯

		ha: 'left'		ha: 'centre'	ha: 'right'
		va: 'centre'	va: 'centre'	va: 'centre'

		 ▯ ▯ ▯							   ▯ ▯ ▯
		 ▯ ▯ ▯ ▯		   Not			 ▯ ▯ ▯ ▯
		 ▯ ▯ ▯ ▯		  Valid			 ▯ ▯ ▯ ▯
		 ▯ ▯ ▯							   ▯ ▯ ▯

		ha: 'left'		ha: 'centre'	ha: 'right'
		va: 'bottom'	va: 'bottom'	va: 'bottom'

		 ▯ ▯ ▯			  ▯ ▯ ▯			   ▯ ▯ ▯
		 ▯ ▯ ▯			  ▯ ▯ ▯			   ▯ ▯ ▯
		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯
		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯		 ▯ ▯ ▯ ▯

	sharex, sharey : bool or {'none', 'all', 'row', 'col'}, default: False
		Not implemented.

	wspace : float, optional
		The horzontal spacing between figure subplots, expressed as a fraction of the subplot width.

	hspace : float, optional
		The vertical spacing between figure subplots, expressed as a fraction of the subplot height.

	sharex, sharey : bool or {'none', 'all', 'row', 'col'}, default: False
		As per matplotlib's usage, controls sharing of properties among x (`sharex`) or y (`sharey`)
		axes:

			- True or 'all': x-/y-axis will be shared among all subplots.
			- False or 'none': each subplot x-/y-axis will be independent (default).
			- 'row': each subplot row will share an x-/y-axis
			- 'col': each subplot column will share an x-/y-axis

		When subplots have a shared x-axis along a column, only the x tick
		labels of the last complete row of the subplot are created. Similarly, 
		when subplots have a shared y-axis along a row, only the y tick labels of the
		first complete column subplot are created. To later turn other subplots'
		ticklabels on, use `~matplotlib.axes.Axes.tick_params`.

	squeeze : bool, optional, default: True
		As per matplotlib's usage, the following applies:
			- If True, extra dimensions are squeezed out from the returned
			  array of Axes:

				- if only one subplot is constructed (nrows=ncols=1), the
				  resulting single Axes object is returned as a scalar.
				- for Nx1 or 1xM subplots, the returned object is a 1D numpy
				  object array of Axes objects.
				- for NxM, subplots with N>1 and M>1 are returned
				  as a 2D array.

			- If False, no squeezing at all is done: the returned Axes object
			  is always a 2D array containing Axes instances, even if it ends
			  up being 1x1.

		If naxes < ncols*nrows, the only sensible option is to return a 1D numpy array
	
	figsize : 2-tuple of floats, default: rcParams["figure.figsize"] * (ncols, nrows)
		The dimensions of the figure (width, height) in inches. If not specified, the default is to scale
		the default rcParams figure.figsize by the number of rows or columns.

	axes_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib `subplot` function.

	**kwargs : Text instance properties
		kwargs are used to specify properties of `subplots` instances
		A list of valid `axis` kwargs can be found here:
		[https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes](https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes "Matplotlib.axes.Axes")
	
	Returns
	-------
	fig : pyplot.Figure
	
	axes : pyplot.axes.Axes object or array of Axes objects
		axes may either be a single Axes object or a numpy array if naxes > 1.
		The dimensions of this array are controlled by the squeeze keyword above.

	"""

	from .base_func import dict_splicer
	
	import warnings
	from matplotlib import rcParams
	from matplotlib.gridspec import GridSpec
	from matplotlib.pyplot import figure, subplot
	
	from numpy import ceil, array, reshape, empty

	gridRef = [[0,0], [1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [2,5], [3,4], [3,4], 
			   [4,4], [3,5], [3,5], [4,4], [3,6], [3,6], [4,5], [4,5], [3,7], [5,5], [5,5], [4,6], [5,5]]


	if (isinstance(sharex, bool)):
		sharex = "all" if sharex else "none"
	if (isinstance(sharey, bool)):
		sharey = "all" if sharey else "none"

	if (naxes == None): # No number of axes specified
		if (nrows==None): nrows=1
		if (ncols==None): ncols=1
		
		naxes = nrows*ncols
	else:
		if (ncols == None and nrows == None):
			if (naxes <= 25):
				nrows, ncols = gridRef[naxes] # Get best combination of rows x cols for number of axes
			else:
				raise NotImplementedError(f"The naxes parameter is currently not implemented for naxes > 25.")
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
	
	# Specify the row/col origin
	row0 = (nrows-1)*2 if va=='bottom' else 0
	col0 = (ncols-1)*2 if ha=='right' else 0

	axes = empty(naxes, dtype=object) # Create the empty array to hold each axis
	for ii in range(naxes):
		row = 2*(nrows-(ii//ncols)-1) if va=='bottom' else 2*(ii//ncols) # current row
		col = 2*(ncols-(ii%ncols)-1) if ha=='right' else 2*(ii%ncols) # current column

		# Select which axes this axis needs to be shared with
		sharewith = {"none": None, "all": axes[0],
					 "row": axes[(ii//ncols)*ncols], "col": axes[ii%ncols]}
		
		axes_kw[ii]["sharex"] = sharewith[sharex]
		axes_kw[ii]["sharey"] = sharewith[sharey]
			
		if (row == (0 if va=='bottom' else (nrows-1)*2) and ha=='centre'):
			axes[ii] = subplot(gs[row:row+2,col+delta:col+delta+2],**axes_kw[ii])
		elif (col == (0 if ha=='right' else (ncols-1)*2) and va=='centre'):
			axes[ii] = subplot(gs[row+delta:row+delta+2,col:col+2],**axes_kw[ii])
		else:
			axes[ii] = subplot(gs[row:row+2,col:col+2],**axes_kw[ii])


	# turn off redundant tick labeling
	if sharex in ["col", "all"]:
		if (ha=='centre'):
			warnings.warn("Removing redundant shared xtick labels not possible when ha='centre'")
		else:
			# turn off all but the bottom row
			for ax in axes[ncols:] if va=='bottom' else axes[:naxes-ncols]:
				ax.xaxis.set_tick_params(which='both',
										 labelbottom=False, labeltop=False)
				ax.xaxis.offsetText.set_visible(False)

	if sharey in ["row", "all"]:
		if (va=='centre'):
			warnings.warn("Removing redundant shared ytick labels not possible when va='centre'")
		else:
			# turn off all but the leftmost column
			for ii, ax in enumerate(axes):
				if (ha=='left'):
					if (ii%ncols!=0):
						ax.yaxis.set_tick_params(which='both',
												 labelbottom=False, labeltop=False)
						ax.yaxis.offsetText.set_visible(False)
				elif (ha=='right'):
					if (ii%ncols+1!=ncols and ii!=naxes-1):
						ax.yaxis.set_tick_params(which='both',
												 labelbottom=False, labeltop=False)
						ax.yaxis.offsetText.set_visible(False)


	### Squeeze axes object
	if (squeeze == True):
		# Discarding unneeded dimensions that equal 1.  If we only have one
		# subplot, just return it instead of a 1-element array.
		return (fig, axes.item()) if axes.size == 1 else (fig, axes.squeeze())
	else:
		# Returned axis array will be always 2-d, even if nrows=ncols=1.
		if (naxes != ncols*nrows):
			warnings.warn("squeeze = False not possible when naxes < nrows*ncols.")
			return (fig, axes.squeeze())
		else:
			return (fig, reshape(axes, (nrows, ncols)))