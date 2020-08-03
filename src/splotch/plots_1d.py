########################################################################
############## Definition of all wrappers for 1D plotting ##############
########################################################################

####################################
# Generalized lines
####################################
def axline(x=None,y=None,a=None,b=None,label=None,lab_loc=0,ax=None,plot_kw={},**kwargs):
	
	"""Generalised axis lines.
	
	This function aims to generalise the usage of axis lines calls (axvline/axhline) together and
	to allow lines to be specified by a slope/intercept according to the function y=a*x + b.
	
	Parameters
	----------
	x : int or list, optional 
		x position(s) in data coordinates for a vertical line(s).
	y : int or list, optional
		y position(s) in data coordinates for a horizontal line(s).
	a : int or list, optional
		Slope(s) of diagonal axis line(s), defaults to 1 if not specified when b is given.
	b : int or list, optional
		Intercept points(s) of diagonal axis line(s), defaults to 0 if not specified when a is given.
	label : str, optional
		Sets label(s) for line(s) and plots legend.
	lab_loc : int, optional
		Defines the position of the legend. Defaults as lab_loc=0.
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
	**kwargs: Line2D properties, optional
		kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
		antialiasing, etc. A list of available `Line2D` properties can be found here: 
		https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
	
	Returns
	-------
	lines
		A list of Line2D objects representing the plotted data.
	
	"""
	
	from matplotlib.pyplot import plot, legend, gca
	from .base_func import axes_handler,plot_finalizer,dict_splicer,is_numeric
	from warnings import warn
	
	# Handle deprecated variables
	deprecated = {'plabel':'label'}
	for dep in deprecated:
		if dep in kwargs:
			warn(f"'{dep}' will be deprecated in future verions, using '{deprecated[dep]}' instead")
			if (dep=='plabel'): label = kwargs.pop(dep)
	
	if ax is not None:
		old_axes=axes_handler(ax)
	else:
		ax=gca()
		old_axes=ax
	
	if not (any([is_numeric(var) for var in [x,y,a,b]])): # If nothing has been specified
		raise TypeError("axline() missing one of optional arguments: 'x', 'y', 'a' or 'b'")
	
	for i, val in enumerate([x,y,a,b]):
		if (val is not None):
			try: # Test whether the parameter is iterable
				temp=(k for k in val)
			except TypeError: # If not, convert to a list
				if   (i == 0): x=[x]
				elif (i == 1): y=[y]
				elif (i == 2): a=[a]
				elif (i == 3): b=[b]
	
	if (x is not None and y is not None): # Check whether both x and y were specified
		raise ValueError("'x' and 'y' cannot be both specified")
	
	if (x is not None): # Check conditions if x specified
		if (any([a,b])): # Should not specify a or b, if x given.
			raise ValueError("'{0}' cannot be specified if x specified".format('a' if a else 'b'))
		L=len(x)
	
	if (y is not None): # Check conditions if y specified
		if (any([a,b])): # Should not specify a or b, if y given.
			raise ValueError("'{0}' cannot be specified if y specified".format('a' if a else 'b'))
		L=len(y)
	
	if (a is not None):
		if (b is None): # If no intercept specified
			b=[0]*len(a) # set b to 0 for all a
		else:
			if (len(b) == 1):
				b=[b[0]]*len(a)
			elif (len(b) != len(a)):
				if (len(a) == 1):
					a=[a[0]]*len(b)
				else:
					raise ValueError(f"Length of 'a' ({len(a)}) and length of 'b' ({len(b)}) must be equal or otherwise 1")
		L=len(a)
	elif (b is not None):
		if (a is None): # If no slope specified
			a=[1]*len(b) # set a to 1 for all b
		L=len(b)
	
	if type(label) is not list:
		label=[label for i in range(L)]
	elif (len(label) != L):
		raise ValueError("Length of label list ({0}) must match the number of lines given ({1}).".format(len(label),L))
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par={**plot_kw, **kwargs} # For Python > 3.5
	plot_par=plot_kw.copy()
	plot_par.update(kwargs)
	
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par=dict_splicer(plot_par,L,[1]*L)
	
	lines=[] # Initialising list which contains each line
	if (x is not None):
		for ii, xx in enumerate(x):
			l=ax.axvline(x=xx,**plot_par[ii],label=label[ii])
			lines.append(l)
	if (y is not None):
		for ii, yy in enumerate(y):
			l=ax.axhline(y=yy,**plot_par[ii],label=label[ii])
			lines.append(l)
	if (a is not None):
		for ii, pars in enumerate(zip(a,b)):
			aa=pars[0]; bb=pars[1]
			
			xLims=ax.get_xlim()
			yLims=ax.get_ylim()
			
			lines.append(plot([xLims[0],xLims[1]],[aa*xLims[0]+bb,aa*xLims[1]+bb],label=label[ii],**plot_par[ii]))
			
			ax.set_xlim(xLims)
			ax.set_ylim(yLims)
			
	if any(label):
		legend(loc=lab_loc)
	if ax is not None:
		old_axes=axes_handler(old_axes)

	return lines[0] if len(lines) == 1 else lines

####################################
# Broken axis plot
####################################
def brokenplot(x,y=None,xbreak=None,ybreak=None,xlim=None,ylim=None,sep=0.05,
			   xinvert=False,yinvert=False,xlog=False,ylog=False,title=None,
			   xlabel=None,ylabel=None,label=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):

	"""Broken Axis Plot Function
	
	Creates a standard plot call with an axis break at `xbreak` or `ybreak` for vertical or
	horizontal breaks.
	
	Parameters
	----------
	x : array-like or list
		If list it is assumed that each elemement is array-like. If y is not given, the given values
		pass to y and a numpy array is generated with numpy.arange() for the x values.
	y : array-like or list, optional
		If list it is assumed that each elemement is array-like.
	xbreak/ybreak : float or tuple-like, required
		The location(s) of the vertical or horizontal breaks is controlled by xbreak or ybreak, respectively.
		The value can be a single location or a tuple defining the (start, stop) points of the break.
		Only one coordinate can be broken in a given plot.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	sep : float, optional, default: 0.05
		The separation size of the axis break, given as a fraction of the axis dimensions.
	xinvert : bool or list, optional
		If true inverts the x-axis.
	yinvert : bool or list, optional
		If true inverts the y-axis.
	xlog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	label : str, optional
		Sets the legend for the plot.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
	**kwargs: Line2D properties, optional
		kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
		antialiasing, etc. A list of available `Line2D` properties can be found here: 
		https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
	
	Returns
	-------
	lines
		A list of Line2D objects (paired as tuples) representing the plotted data.
		The lines are given as pairs to correspond to the separate lines either side of the x/ybreak.
	
	"""
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	from numpy import shape, arange, ndarray
	from matplotlib.pyplot import plot, legend, show, sca, gca
	from matplotlib.transforms import Bbox
	from warnings import warn
	
	# Handle deprecated variables
	deprecated = {'plabel':'label'}
	for dep in deprecated:
		if dep in kwargs:
			warn(f"'{dep}' will be deprecated in future verions, using '{deprecated[dep]}' instead")
			if (dep=='plabel'): label = kwargs.pop(dep)
	
	if ax is not None:
		old_axes=axes_handler(ax)
	else:
		ax=gca()
		old_axes=ax

	if type(x) is not list or len(shape(x))==1:
		x=[x]
	L=len(x)
	if y is None:
		y=x
		x=[arange(len(x[i])) for i in range(L)]
	else:
		if type(y) is not list or len(shape(y))==1:
			y=[y]
	if type(label) is not list:
		label=[label for i in range(L)]
	
	# Validate x/ybreak
	if (xbreak == None):
		raise ValueError("Require either xbreak/ybreak to be specified.")

	if (ybreak != None):
		raise NotImplementedError("ybreak not yet implemented.")


	if type(xbreak) not in [list,tuple,ndarray]:
		xbreak=(xbreak, xbreak)
	else:
		if (len(xbreak) != 2):
			raise ValueError("xbreak must be a single value of a tuple-like list of two elements.")
	
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par={**plot_kw, **kwargs} # For Python > 3.5
	plot_par=plot_kw.copy()
	plot_par.update(kwargs)
	
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par=dict_splicer(plot_par,L,[1]*L)
	
	# Get the original axis position
	pos0=ax.get_position(original=True)
	width0, height0=pos0.x1 - pos0.x0, pos0.y1 - pos0.y0
	
	lines=[] # Initialising list which contains each line
	for i in range(L):
		# First side plot call
		l1=ax.plot(x[i],y[i],label=label[i],**plot_par[i])
		
		# Get the axis limits if not already specified
		xlims=ax.get_xlim() if xlim == None else xlim
		ylims=ax.get_ylim() if ylim == None else ylim
	
		# Define the positions of the two separated axes
		if (i == 0):
			pos1=Bbox(list(pos0.get_points()))
			pos1.x1=pos1.x0 + (pos1.x1-pos1.x0)*(sum(xbreak)/2-xlims[0])/(xlims[1]-xlims[0]) - sep*(pos1.x1-pos1.x0)/2
			
			pos2=Bbox(list(pos0.get_points()))
			pos2.x0=pos2.x0 + (pos2.x1-pos2.x0)*(sum(xbreak)/2-xlims[0])/(xlims[1]-xlims[0]) + sep*(pos2.x1-pos2.x0)/2
			
			ax.set_position(pos1) # Resize the first axis
			ax2=ax.figure.add_axes(pos2) # Add and duplicate the plotting in the second axis
			
			# Set the new axis limits at the break point
			ax.set_xlim(xlims[0],xbreak[0])
			ax2.set_xlim(xbreak[1],xlims[1])
		
		# Second side plot call
		l2=ax2.plot(x[i],y[i],label=None,**plot_par[i])

		lines.append((*l1,*l2)) # Add line as tuple of both axes.
		
		width1, height1=pos1.x1 - pos1.x0, pos1.y1 - pos1.y0
		width2, height2=pos2.x1 - pos2.x0, pos2.y1 - pos2.y0
		
		dx1, dy1=0.01 * width0/(width0-width1-sep/2), height1*0.025
		
		dash_kw=dict(transform=ax2.transAxes, color='black', linestyle='-', marker='', clip_on=False)
		ax2.plot((0 - dx1, 0 + dx1), (0 - dy1, 0 + dy1), **dash_kw)  # bottom-right diagonal
		ax2.plot((0 - dx1, 0 + dx1), (1 - dy1, 1 + dy1), **dash_kw)  # top-right diagonal
		
		dx2, dy2=0.01 * width0/(width0-width2-sep/2), height2*0.025
		dash_kw.update(transform=ax.transAxes)  # switch to the left axes
		ax.plot((1 - dx2, 1 + dx2), (0 - dy2, 0 + dy2), **dash_kw)  # bottom-left sep/5iagonal
		ax.plot((1 - dx2, 1 + dx2), (1 - dy2, 1 + dy2), **dash_kw)  # top-left sep/5iagonal
		
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelright=False,which='both')  # don't put tick labels at the top
		ax.yaxis.tick_left()
		
		ax2.spines['left'].set_visible(False)
		ax2.tick_params(labelleft=False,which='both')  # don't put tick labels at the top
		ax2.yaxis.tick_right()

		# Check that there is no duplicate ticks over both axes	
	if (xbreak):
			if (ax.get_xticks()[-1] == ax2.get_xticks()[0]):
				if (xbreak[0] >= (xlims[0] + xlims[1])*0.5):
					ax.set_xticks(ax.get_xticks()[:-1]) # Remove duplicate tick on left side
				else:
					ax2.set_xticks(ax2.get_xticks()[1:]) # Remove duplicate tick on right side
	sca(ax)
	
	if any(label):
		ax.legend(loc=lab_loc)
	
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	
	if ax is not None:
		old_axes=axes_handler(old_axes)

	return (lines[0] if len(lines) == 1 else lines)

####################################
# Curves from analytic expressions
####################################
def curve(expr, var=None, subs={}, permute=False, bounds=None, num=101, xlim=None, ylim=None, xinvert=False, yinvert=False,
		  xlog=False, ylog=False, title=None, xlabel=None, ylabel=None, label=True, lab_loc=0,
		  grid=None, ax=None, plot_kw={}, **kwargs):
	"""Function Plotting
	
	Plot the curve corresponding to a definingned function over the range of [from, to].
	
	Parameters
	----------
	expr : str, sympy.Expr or callable()
		An expression parsed either as a string, sympy expression or callable (function or lambda)
		which will be evaluated by the function in the range of `bounds`.
	var : str or sympy symbol, default: 'x'
		The independent variable on which to evaluate the expression (i.e. the variable on the x-axis).
		This defaults to the first non-numeric element of the expression or otherwise simply assumes
		this to be 'x'.
	subs : dict, optional
		If `expr` contains more symbols than the independent variable `var`, this dictionary will
		substitute numerical values for all additonal symbols given. `subs` is required if
		additional symbols are specified.
	bounds : float, optional
		The range over which the function will be plotted. If not given, these default to
		the current bounds of the plot.
	num : int, optional (default: 101)
		The number of values along the independent variable on which to evaulate `expr`.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	xinvert : bool or list, optional
		If true inverts the x-axis.
	yinvert : bool or list, optional
		If true inverts the y-axis.
	xlog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	label : bool, str or list-like, optional (default: True)
		Sets the label(s) of the curve(s) for the legend. `label` can be one of:
			- True:
				Creates a label for every curve defined by `subs`. The label will list all values
				`subs` that produced the curve. If a parameter in `subs` has only one value
				(i.e. constant amongst all curves), it will not appear in the label.
			- str:
				If a single string is given, only one label will be shown and all curves will
				be shown as the label handle.
			- list (length must equal number of curves):
				A label given to each curve that is produced.
			- False/None:
				No label will be assigned and a legend will not be shown.
	lab_loc : int, optional
		Defines the position of the legend
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D
		properties. It is recommended that kwargs be parsed implicitly through **kwargs
		for readability.
	**kwargs: Line2D properties, optional
		kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, 
		antialiasing, etc. A list of available `Line2D` properties can be found here: 
		https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
	
	Returns
	-------
	curves : list of (or single) pyplot.Line2D object(s)
		A list of Line2D objects created for each curved create by `subs`.
	expr : Sympy.Expr
		If expr was given as a string, this returns the sympy expression created from `sympy.sympify()`.
		Otherwise, simply returns the `expr` that was given.
	
	"""
	
	from .base_func import axes_handler,dict_splicer,plot_finalizer,simpler_dict_splicer
	
	from sympy import symbols, sympify, Expr
	from sympy.utilities.lambdify import lambdify
	from numpy import linspace, logspace, log10, empty, array, meshgrid, prod
	from collections import Iterable
	
	from matplotlib.pyplot import plot, legend, gca
	from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D, HandlerTuple
	from matplotlib import rcParams
	from warnings import warn
	
	# Handle deprecated variables
	deprecated = {'plabel':'label'}
	for dep in deprecated:
		if dep in kwargs:
			warn(f"'{dep}' will be deprecated in future verions, using '{deprecated[dep]}' instead")
			if (dep=='plabel'): label = kwargs.pop(dep)
	
	if ax is not None:
		old_axes=axes_handler(ax)
	else:
		ax=gca()
		old_axes=ax
	
	isfunc=False
	if (isinstance(expr, str)):
		expr=sympify(expr)
	elif (callable(expr)):
		isfunc=True
	elif (isinstance(expr, Expr)):
		pass
	else:
		raise TypeError(f"`expr` must be of type `str` or sympy.Expr, instead got {type(expr)}.")
	
	if (isfunc == False):
		symbs=expr.free_symbols # Get the symbols in the expression
		symbkeys=[str(sym) for sym in symbs]
		
		if (var != None): # Validate the independent variable
			if (var not in symbkeys):
				raise ValueError(f"Independent variable '{var}' was not found in 'expr'.")
		else: # Assume independent variable is 'x', otherwise, assume the first symbol.
			var='x' if 'x' in symbkeys or len(symbkeys)==0 else symbkeys[0]
		
		# Validate the substitution variable names
		if (subs != None):
			if (var in list(subs)):
				raise ValueError(f"Independent variable '{var}' should not be in subs.")
			for key in list(subs):
				if (key not in symbkeys):
					raise KeyError(f"Substitution variable '{key}' does not exist in 'expr'.")
	
	# The lengths of each substitute value list, len=1 if just a single value
	lens=[len(subs[key]) if (isinstance(subs[key], Iterable) and type(subs[key])!=str) else 1 for key in list(subs)]
	
	if (permute == True):
		L=prod(lens)
		perms=array(meshgrid(*subs.values())).reshape(len(subs),-1)
		permsubs={}
		for ii, key in enumerate(list(subs)):
			permsubs[key]=perms[ii]
		subsarr=simpler_dict_splicer(permsubs,L,[1]*L)
	else:
		L=max(lens) if len(lens) > 0 else 1
		subsarr=simpler_dict_splicer(subs,L,[1]*L)
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par={**plot_kw, **kwargs} # For Python > 3.5
	plot_par=plot_kw.copy()
	plot_par.update(kwargs)
	
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par=dict_splicer(plot_par,L,[1]*L)
	
	if (bounds == None):
		if (xlim != None):
			bounds=xlim
		else:
			bounds=ax.get_xlim()
	
	vararr=logspace(*log10(bounds),num=num) if xlog else linspace(*bounds,num=num)
	
	curves=[None]*L
	for ii in range(L):
		if (isfunc):
			curvearr=expr(vararr, **subsarr[ii])
		else:
			func=lambdify(var, expr.subs(subsarr[ii]), 'numpy') # returns a numpy-ready function
			curvearr=func(vararr)
		curves[ii]=plot(vararr,curvearr,**plot_par[ii])[0]
	
	# Create the legend object
	if (label == True):
		labellist=['']*L
		for ii in range(L): # Make a label for each of sub values
			# Create a list of the sub names and their values.
			labellist[ii]="; ".join([f"{key}={subsarr[ii][key]}" for jj, key in enumerate(list(subsarr[ii]))]) # join substitute strings together
		# else:
		# 	print([f"{key}={subsarr[0][key]}" for jj, key in enumerate(list(subsarr[0]))])
		# 	labellist="; ".join([f"{key}={subsarr[0][key]}" for jj, key in enumerate(list(subsarr[0]))]) # join substitute strings together
		
		ax.legend(handles=curves, labels=labellist, loc=lab_loc)
	
	elif (isinstance(label, Iterable)): # A list-like iterable object or string has been given
		if (type(label) == str): # Single string given, multiple lines for the label.
			ax.legend(handles=[tuple(curves)], labels=[label], loc=lab_loc,handler_map={tuple: HandlerTuple(ndivide=L,pad=1/L**0.8)},
						handlelength=rcParams["legend.handlelength"]*L**0.8)
		else:
			if (len(label) == L): # number of labels matches number of curves
				labellist=label
			else:
				raise ValueError("Length of 'label' list does not match the number of curves defined by 'subs'.")
			ax.legend(handles=curves, labels=labellist, loc=lab_loc)
	elif (label == False or label == None): # Do not plot a legend
		pass
	else:
		raise ValueError("Value of 'label' is invalid. Expected bool, str, iterable or None, but got '{0}'.".format(label))
	
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	
	if ax is not None:
		old_axes=axes_handler(old_axes)
	
	return(curves[0] if len(curves)==1 else curves, expr)

####################################
# 1D histogram and binned statistics
####################################
def hist(data,bin_type=None,bins=None,dens=True,cumul=None,scale=None,weights=None,hist_type=None,v=None,vstat=None,
			xlim=None,ylim=None,nmin=0,xinvert=False,yinvert=False,xlog=False,ylog=None,title=None,xlabel=None,ylabel=None,
			label=None,lab_loc=0,ax=None,grid=None,plot_kw={},output=None,**kwargs):
	
	"""1D histogram function.
	
	The plotting is done with pyplot.plot(), so histograms are shown with interpolated curves instead
	of the more common stepwise curve. For this reason splotch.histstep is a better choice for
	small datasets. 
	
	Parameters
	----------
	data : array-like or list
		If list it is assumed that each elemement is array-like.
	bin_type : {'number','width','edges','equal'}, optional
		Defines how is understood the value given in bins: 'number' for the desired number of bins,
		'width' for the width of the bins, 'edges' for the edges of bins, and 'equal' for making
		bins with equal number of elements (or as close as possible). If not given it is inferred
		from the data type of bins: 'number' if int, 'width' if float and 'edges'if ndarray.
	bins : int, float, array-like or list, optional
		Gives the values for the bins, according to bin_type.
	dens :  bool or list, optional
		If false the histogram returns raw counts.
	cumul : bool or list, optional
		If true, produces a cumulative distribution instead of a histogram.
	scale : float or list, optional
		Scaling to be applied to the counts.
	weights : array-like or None, optional
		An array of weights with the same shape as data. For each value in data, it will only
		contribute its given weight towards the bin count (instead of 1). If dens is True, weights
		will also be normalised so that the integral over the density remains 1. Default: None
	hist_type : str, optional.
		Defines the type of histogram to be drawn. 'line' and 'step' produce lines, with the former
		drawing lines conecting the values of each bin positioned on their centre, and the latter
		drawing a stepwise line, with the edges of each step coinciding with the bin edges.
		'bar' produces a bar plot. All have filled version (i.e., 'linefilled'), which fills the
		space between the edges of the histogram and 0.
	v : array-like or list, optional
		If a valid argument is given in vstat, defines the value used for the binned statistics.
	vstat : str, function  or list, optional
		Must be or contain one of the valid str arguments for the statistics variable
		in scipy.stats.binned_statistic ('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or
		function(s) that takes a 1D array and outputs an integer or float.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	nmin : int, optional (default: 0)
		The minimum number of points required in a bin in order to be plotted.
	xinvert : bool or list, optional
		If true inverts the x-axis.
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	yinvert : bool or list, optional
		If true inverts the y-axis.
	xlog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	label : str, optional
		Sets the label for the plot.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	plot_par : dict, optional
		Passes the given dictionary as a kwarg to the plotting function.
	output : boolean, optional
		If True, returns the edges and values of the histogram.
	
	Returns
	-------
	n : list
		List containing the arrays with the values for each histogram drawn. Only provided
		if output is True.
	bins_edges : list
		List containing the arrays with the bin edges for each of the histograms drawn.
		Only provided if output is True.
	"""
	
	from numpy import cumsum as np_cumsum, sum as np_sum, shape
	from scipy.stats import binned_statistic
	from numpy import array, ndarray, diff, dtype, histogram, inf, nan, nanmax, nanmean, nanstd, ones, where
	from matplotlib.pyplot import bar, fill_between, gca, legend, plot, rcParams, step
	from .base_func import axes_handler,bin_axis,dict_splicer,plot_finalizer,step_filler
	from warnings import warn
	
	# Handle deprecated variables
	deprecated = {'plabel':'label'}
	for dep in deprecated:
		if dep in kwargs:
			warn(f"'{dep}' will be deprecated in future verions, using '{deprecated[dep]}' instead")
			if (dep=='plabel'): label = kwargs.pop(dep)

	if ax is not None:
		old_axes=axes_handler(ax)
	if type(data) not in [list, tuple, ndarray] or (len(shape(data)) == 1 and array(data).dtype is not dtype('O')):
		data=[data]
	L=len(data)
	if type(bin_type) not in [list, tuple]:
		bin_type=[bin_type]*L
	if type(bins) not in [list, tuple, ndarray]:
		if bins is not None:
			bins=[bins]*L
		else:
			bins=[int((len(d))**0.4) for d in data]
	if type(weights) not in [list, tuple, ndarray] or (len(shape(weights)) == 1):
		weights=[weights]*L
	if type(dens) not in [list, tuple]:
		dens=[dens]*L
	if type(cumul) not in [list, tuple]:
		cumul=[cumul]*L
	if type(scale) not in [list, tuple, ndarray]:
		scale=[scale]*L
	if type(v) not in [list, tuple, ndarray] or (len(shape(v)) == 1):
		v=[v]*L
	if type(nmin) not in [list, tuple, ndarray] or (len(shape(nmin)) == 1):
		nmin=[nmin]*L
	if type(vstat) not in [list, tuple]:
		vstat=[vstat]*L
	if type(label) not in [list, tuple]:
		label=[label for i in range(L)]
	if type(xlim) in [int,float]:
		xlim=[nanmean(data)-xlim*nanstd(data),nanmean(data)+xlim*nanstd(data)]
	
	if None in [ylog,hist_type,output]:
		from .defaults import Params
		if ylog is None:
			ylog=Params.hist1D_yaxis_log
		if hist_type is None:
			hist_type=Params.hist1D_histtype
		if output is None:
			output=Params.hist1D_output
	if type(hist_type) not in [list, tuple, ndarray]:
		hist_type=[hist_type]*L
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par={**plot_kw, **kwargs} # For Python > 3.5
	plot_par=plot_kw.copy()
	plot_par.update(kwargs)
	# Check if width is given as a kwarg
	if 'width' in plot_par.keys():
		import warnings
		warnings.warn('Received kwarg width, this will be ignored in the histogram',UserWarning)
		if hist_type!='bar':
			temp=plot_par.pop('width')
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par=dict_splicer(plot_par,L,[1]*L)
	
	plot_type={'line':plot,'linefilled':fill_between,'step':step,'stepfilled':step_filler,'bar':bar,'barfilled':bar}
	hist_centre={'line':True,'linefilled':True,'step':False,'stepfilled':False,'bar':False,'barfilled':False}
	bin_edges=[]
	n_return=[]
	
	for i in range(L):
		temp_data,bins_hist,bins_plot=bin_axis(data[i],bin_type[i],bins[i],log=xlog,plot_centre=hist_centre[hist_type[i]])
		n_check=histogram(temp_data,bins=bins_hist,density=False)[0]
		n_check=n_check>=nmin[i]
		if vstat[i]:
			temp_y=binned_statistic(temp_data,v[i],statistic=vstat[i],bins=bins_hist)[0]
		else:
			temp_y=histogram(temp_data,bins=bins_hist,density=dens[i],weights=weights[i])[0]
		if cumul[i]:
			temp_y=np_cumsum(temp_y)
			if dens[i]:
				temp_y=temp_y.astype('float')/nanmax(temp_y)
		if scale[i]:
			if dens[i]:
				temp_y*=len(data[i])/scale[i]
			else:
				temp_y=temp_y.astype('float')/scale[i]
		if ylog:
			temp_y=where(temp_y==0,nan,temp_y)
		temp_y=where(n_check,temp_y,nan)
		y=temp_y
		if hist_type[i]=='step':
			if ylog:
				y=array([y[0]]+[j for j in y])
			else:
				bins_plot=array([bins_plot[0]]+[b for b in bins_plot]+[bins_plot[-1]])
				y=array([0,y[0]]+[j for j in y]+[0])
		if 'bar' in hist_type[i]:
			#prop_cycle=rcParams['axes.prop_cycle']
			#barcolor=prop_cycle.by_key()['color']
			plot_par[i]['width']=diff(bins_plot)
			bins_plot=(bins_plot[0:-1]+bins_plot[1:])/2
			if hist_type[i]=='bar':
				if 'edgecolor' not in plot_par[i].keys():
					p=plot(bins_plot[0],0)
					plot_par[i]['edgecolor']=p[0].get_color()
					p.pop()
					temp_ax=gca()
					temp_ax.relim()
					temp_ax.autoscale()
				plot_par[i]['fill']=False
		plot_type[hist_type[i]](bins_plot,y,label=label[i],**plot_par[i])
		bin_edges.append(bins_plot)
		n_return.append(temp_y)
	if any(label):
		legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)
	if len(n_return)==1:
		n_return=n_return[0]
	if len(bin_edges)==1:
		bin_edges=bin_edges[0]
	if output:
		return(n_return,bin_edges)

####################################
# Standard plots
####################################
def plot(x,y=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,title=None,xlabel=None,
			ylabel=None,label=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""Base plotting function.
	
	This is a wrapper for pyplot.plot().
	
	Parameters
	----------
	x : array-like or list
		If list it is assumed that each elemement is array-like. If y is not given, the given
		values pass to y and a numpy array is generated with numpy.arange() for the x values.
	y : array-like or list, optional
		If list it is assumed that each elemement is array-like.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	xinvert : bool or list, optional
		If true inverts the x-axis.
	yinvert : bool or list, optional
		If true inverts the y-axis.
	xlog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool or list, optional
		If True the scale of the x-axis is logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	label : str, optional
		Sets the legend for the plot.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are
		Line2D properties.
	**kwargs: Line2D properties, optional
		kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
		antialiasing, etc. A list of available `Line2D` properties can be found here: 
		https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
	
	Returns
	-------
	lines
		A list of Line2D objects representing the plotted data.
	"""
	
	from numpy import shape, arange
	from matplotlib.pyplot import plot, legend
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	from warnings import warn

	# Handle deprecated variables
	deprecated = {'plabel':'label'}
	for dep in deprecated:
		if dep in kwargs:
			warn(f"'{dep}' will be deprecated in future verions, using '{deprecated[dep]}' instead")
			if (dep=='plabel'): label = kwargs.pop(dep)
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list or len(shape(x))==1:
		x=[x]
	L=len(x)
	if y is None:
		y=x
		x=[arange(len(x[i])) for i in range(L)]
	else:
		if type(y) is not list or len(shape(y))==1:
			y=[y]
	if type(label) is not list:
		label=[label for i in range(L)]
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par={**plot_kw, **kwargs} # For Python > 3.5
	plot_par=plot_kw.copy()
	plot_par.update(kwargs)
	
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par=dict_splicer(plot_par,L,[1]*L)
	
	lines=[] # Initialising list which contains each line
	for i in range(L):
		lines += plot(x[i],y[i],label=label[i],**plot_par[i])
	if any(label):
		legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)
	
	return (lines[0] if len(lines) == 1 else lines)
