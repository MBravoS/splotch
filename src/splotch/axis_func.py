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


def colorbar(mappable=None,ax=None,label='',orientation='vertical',loc=1,transform=None,
			  inset=False,aspect=0.05,width=None,height=None,pad=0.05,ticks=None,bar_kw={},**kwargs):
	"""Colorbar function
	
	This function will produce a colorbar for any currently plotted mappable on a given axis.
	
	Parameters
	----------
	mappable : ScalarMappable, optional
		A list or individual matplotlib.cm.ScalarMappable (i.e., Image, ContourSet, etc.) described by this colorbar(s).
		This argument is optional and will seek out any mappables currently present in each axis given if
		nothing is specified.	
	ax : pyplot.Axes, optional
		Use the given axes to produce the colorbars onto. If multiple axes are given, number of mappable objects given must be 
		either one or equal to the number of axis objects. If no mappables are provided, a mappable will be searched for
		independently for each axis. Defaults to the current axis.
	label : str, optional
		The label to be given to the colorbar
	orientation : optional
		The orientation of the colorbar specified either as 'vertical' | 'horizontal'.
		Orientation is necessary to decide which axis of the colorbar to place labels. Default: 'vertical'.
	loc : int or tuple-like, optional
		Specifies the location of the colorbar. Can be a two element tuple-like object in the format
		of (x0, y0).
	transform : matplotlib.transforms.Transform instance, optional
		The transformation instance to be used for colorbar location if loc is tuple-like. For example, using ax.transAxes()
		will specify the colorbar location in the coordinates of the axis; (0,0) is bottom-left and (1,1)
		is top-right. Default: ax.transAxes.
	inset : boolean, optional
		Whether to inset the colorbar within the inside of the axis. If loc not tuple-like, this will add
		padding to both sides of the colorbar to avoid colliding with the axis spine. Default: False.
	aspect : float, optional
		The aspect ratio of the colorbar always taken as the ratio of the long-side to the short-side.
		Default: 0.05
	width, height : float, optional
		The width and height of the colorbar in the coordinate system of specified `transform`, which by
		default is in the coordinates of the axis.
	pad : float, optional
		The padding given to the colorbar axis offset from the margin of the axis. If loc is specified to an edge
		which has an axis label/ticks, padding will be added from the edge of the label.
	ticks : None, list-like or Locator() object, optional
		Specifies the locations of ticks on th colorbar axis. If None, ticks are determined automatically from the input.
	bar_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are colorbar properties.
	**kwargs : Colorbar properties, optional
		Keyword arguments are used to specify matplotlib.pyplot.colorbar specific properties such as
		extend, spacing, format, drawedges, etc. The list of available properties can be found here: 
		https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.colorbar.html#matplotlib.pyplot.colorbar
	
	Returns
	-------
	cbar : 
		The colorbar object.
	"""
	
	from .base_func import axes_handler, dict_splicer

	from matplotlib.pyplot import gca
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
	from mpl_toolkits.axes_grid1.colorbar import colorbar
	

	# Validate axis input
	if ax is not None:
		try: # check if iterable
			_ = (i for i in ax)
			old_axes = axes_handler(ax[0])
			axes = ax
		except (TypeError):
			old_axes = axes_handler(ax)
			axes = [ax]
			
	else:
		axes = [gca()]
		old_axes = axes[0]

	if type(loc) != int:
		raise NotImplementedError("loc must be specified as integer. Providing loc as a tuple-like as colorbar anchor position is not yet implemented.")

	# Validate aspect value
	if (aspect <= 0 or aspect > 1):
		raise ValueError("Value for aspect must be strictly positive and less than or equal to 1 (i.e. 0 < aspect <= 1)")
	
	# Get width/height values of colorbar
	if (width == None and height == None):
		width, height = (aspect, 1.0-ins*2*pad) if orientation is 'vertical' else (1.0-ins*2*pad, aspect)
	else:
		if height == None:
			height = width/aspect if orientation is 'vertical' else aspect*width
		elif width == None:
			width = aspect*height if orientation is 'vertical' else height/aspect

	### Define the positions of preset colorbars
	labpad = 0.125 # The padding added for labels
	ins = 1 if inset == True else 0 # Convert inset boolean to a binary multiplier
	
	# Vertically-oriented colorbars
	vertPositions = {1:(1+pad-ins*(2*pad+width), 1-height-ins*pad, width, height),
					 2:(0-width-labpad-pad+ins*(labpad+2*pad+width), 1-height-ins*pad, width, height),
					 3:(0-width-labpad-pad+ins*(labpad+2*pad+width), 0+ins*pad, width, height),
					 4:(1+pad-ins*(2*pad+width), 0+ins*pad, width, height),
					 5:(1+pad-ins*(2*pad+width), 0.5*(1-height), width, height),
					 6:(0.5*(1-width), 1+pad-ins*(2*pad+height), width, height),
					 7:(0-width-labpad-pad+ins*(labpad+2*pad+width), 0.5*(1-height), width, height),
					 8:(0.5*(1-width),0-height-labpad-pad+ins*(labpad+2*pad+height), width, height),
					 9:(0.5*(1-width), 0.5*(1-height), width, height)}
	
	# horizontally-oriented colorbars
	horPositions  = {1:(1-width-ins*pad, 1+pad-ins*(2*pad+height), width, height),
					 2:(0+ins*pad, 1+pad-ins*(2*pad+height), width, height),
					 3:(0+ins*pad, 0-height-labpad-pad+ins*(height+labpad+2*pad), width, height),
					 4:(1-width-ins*pad, 0-height-labpad-pad+ins*(height+labpad+2*pad), width, height),
					 5:(1+pad-ins*(2*pad+width), 0.5*(1-height), width, height),
					 6:(0.5*(1-width), 1+pad-ins*(2*pad+height), width, height),
					 7:(0-width-labpad-pad+ins*(width+labpad+2*pad), 0.5*(1-height), width, height),
					 8:(0.5*(1-width), 0-height-labpad-pad+ins*(labpad+2*pad+height), width, height),
					 9:(0.5*(1-width), 0.5*(1-height), width, height)}


	# Combine the `explicit` bar_kw dictionary with the `implicit` **kwargs dictionary
	#bar_par = {**bar_kw, **kwargs} # For Python > 3.5
	bar_par = bar_kw.copy()
	bar_par.update(kwargs)

	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	bar_par = dict_splicer(bar_par,len(axes),[1]*len(axes))

	cbars = [] # Initiate empty list of colorbars for output
	for ii, ax in enumerate(axes):
		if mappable == None:
			mappables = [child for child in ax.get_children() if hasattr(child, 'autoscale_None')]
			if len(mappables) == 1:
				mapper = mappables[0]
			else:
				mapper = mappables[0]
		else:
			try: # check if iterable
				if (len(mappable) == len(axes)):
					mapper = mappable[ii]
				elif (len(mappable) == 1):
					mapper = mappable[0]
				else:
					raise ValueError("Number of mappables given must be either 1 or equal to the number of axes specified.")
			except (TypeError):
				mapper = mappable

		cax = inset_axes(ax, width='100%', height='100%',
						 bbox_to_anchor=vertPositions[loc] if orientation is 'vertical' else horPositions[loc],
						 bbox_transform=transform if transform != None else ax.transAxes,
						 borderpad=0)
		
		cbar = colorbar(mapper, cax=cax, orientation=orientation,ticks=ticks,**bar_par[ii])
		
		# Orient tick axes correctly
		if (orientation is 'horizontal'):
			cbar.ax.set_xlabel(label)
			
			flip = True if loc in [3,4,8] else False # Flip the labels if colorbar on bottom edge
			if (inset==True) and loc not in [5,7,9]: flip = not flip # Reverse flipping in the case of a inset colorbar	
			
			if (flip == True):
				cbar.ax.xaxis.set_label_position('bottom')
				cbar.ax.xaxis.tick_bottom()
			else:
				cbar.ax.xaxis.set_label_position('top')
				cbar.ax.xaxis.tick_top()
		else:
			cbar.ax.set_ylabel(label)
			
			flip = True if loc in [2,3,7] else False  # Flip the labels if colorbar on left edge
			if (inset==True) and loc not in [6,8,9]: flip = not flip # Reverse flipping in the case of a inset colorbar
			
			if (flip == True):
				cbar.ax.yaxis.set_label_position('left')
				cbar.ax.yaxis.tick_left()

		cbars.append(cbar)

	return (cbars[0] if len(cbars) == 1 else cbars)

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