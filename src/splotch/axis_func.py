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

	import matplotlib.text.Text as Text
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if (ax == None):
		ax = plt.gca()

	# Validate `which` value(s)
	whichRef = ['x','y','t','s','l','c',
				'xlabel','ylabel','title','suptitle','legend','colorbar']

	for w in which:
		try:
			wInd = whichRef.index(w)
			wComp = (wInd + len(whichRef)//2) % len(whichRef) # get the complimenting short/long version
			if (whichRef[wComp] in which):
				raise TypeError("adjust_text() received equivalent values for 'which': '{0}' and '{1}'.".format(whichRef[wInd],whichRef[wComp]))

		except (ValueError):
			raise TypeError("adjust_text() received invalid value for 'which' ('{0}'). Must be one of: {1}".format(w,', '.join(whichRef)))
		

	L = len(which)
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	textpar = text_kw.copy()
	textpar.update(kwargs)

	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	textpar = dict_splicer(textpar,L,[1]*L)

	if type(ax) is not list or len(np.shape(ax))==1:
		ax=[ax]

	for a in ax:
		for ii, lab in enumerate(which):
			if (lab in ['x','xlabel']):
				text = a.xaxis.label
			elif (lab in ['y','ylabel']):
				text = a.yaxis.label
			elif (lab in ['t','title']):
				text = a.title
			elif (lab in ['s','suptitle']): # Not implemented
				text = a.title
			elif (lab in ['l','legend']):
				text = a.legend().get_texts() # Get list of text objects in legend
			elif (lab in ['c','colorbar']):
				# Get the axis with the largest ratio between width or height
				caxInd = np.argmax([np.max([c.get_position().width/c.get_position().height,c.get_position().height/c.get_position().width]) for c in axes.figure.get_axes()])
				if (axes.figure.get_axes()[caxInd].get_position().height > axes.figure.get_axes()[caxInd].get_position().width)
					text = ax.figure.get_axes()[caxInd].yaxis.label
				else:
					text = ax.figure.get_axes()[caxInd].xaxis.label

		if (type(text) == Text):
			text.set(**textpar[ii])
		else:
			for t in text:
				t.set(**textpar[ii])
	
	return(None)