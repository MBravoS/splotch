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