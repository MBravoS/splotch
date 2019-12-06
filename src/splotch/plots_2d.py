#### Definition of all wrappers for 2D plotting

#Level contours
def cont(z,x=None,y=None,filled=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,title=None,xlabel=None,
			ylabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""Level contour plotting function (legacy name).
	
	This provides legacy compatibily for old code using the original name of the function.
	This function will eventually be removed, so consider switching to plots_2d.contour().
	"""
	
	import warnings
	
	warnings.warn('This function provides legacy support and it will be removed in the future. Use plots_2d.contour() instead.',FutureWarning)
	
	contour(z,x,y,filled,xlim,ylim,xinvert,yinvert,xlog,ylog,title,xlabel,ylabel,lab_loc,ax,grid,plot_kw={},**kwargs)

def contour(z,x=None,y=None,filled=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,title=None,xlabel=None,
			ylabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""Level contour plotting function.
	
	This is a wrapper for pyplot.contour() and pyplot.contourf().
	
	Parameters
	----------
	z : array-like
		The height values to draw the contours.
	x : array-like, optional
		Position of data points in the x axis.
	y : array-like, optional
		Position of data points in the y axis.
	filled: boolean, optional
		If True draws filled contours. If not given defaults to the value defined in splotch.Params.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	xinvert : bool, optional
		If true inverts the x-axis.
	yinvert : bool, optional
		If true inverts the y-axis.
	xlog : bool, optional
		If True the scale of the x-axis is logarithmic. If not given defaults to the value defined in splotch.Params.
	ylog : bool, optional
		If True the scale of the x-axis is logarithmic. If not given defaults to the value defined in splotch.Params.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	output : boolean, optional
		If True, returns the edges and values of the underlying histogram plus the levels of the contours.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are QuadContourSet properties.
	**kwargs: QuadContourSet properties, optional
		kwargs are used to specify matplotlib specific properties such as cmap, linewidths, hatches, etc.
		The list of available properties can be found here: 
		https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.contour.html
	
	Returns
	-------
	bin_edges_x : array
		The bin edges for the x axis.
	bin_edges_y : array
		The bin edges for the y axis.
	n : array
		The values of the underlying histogram.
	l : array
		The levels for the contours.
	"""
	
	from numpy import shape, linspace
	from matplotlib.pyplot import contour,contourf, legend
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if filled is None:
		from .defaults import Params
		filled=Params.cont_filled
	if x is None:
		x=linspace(0,1,z.shape[0])
	if y is None:
		y=linspace(0,1,z.shape[1])
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)
	
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	
	plotf={False:contour,True:contourf}
	plotf[filled](x,y,z,**plot_par)
	
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

#Errorbars
def errbar(x,y,xerr=None,yerr=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,
	title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""Errorbar plotting function (legacy name).
	
	This provides legacy compatibily for old code using the original name of the function.
	This function will eventually be removed, so consider switching to plots_2d.errorbar().
	"""
	
	import warnings
	
	warnings.warn('This function provides legacy support and it will be removed in the future. Use plots_2d.errorbar() instead.',FutureWarning)
	
	errorbar(x,y,xerr,yerr,xlim,ylim,xinvert,yinvert,xlog,ylog,title,xlabel,ylabel,plabel,lab_loc,ax,grid,plot_kw,**kwargs)

def errorband(x,y,bin_type=None,bins=None,line_stat='mean',band_stat_low='std',band_stat_high='std',line=False,xlim=None,ylim=None,
				xinvert=False,yinvert=False,xlog=False,ylog=None,title=None,xlabel=None,ylabel=None,
				plabel=None,lab_loc=0,ax=None,grid=None,line_kw={},band_kw={},**kwargs):
	
	"""Errorbar plotting function.
	
	This is a wrapper for pyplot.errorbar().
	
	Parameters
	----------
	x : array-like or list
		If list it is assumed that each elemement is array-like.
	y : array-like or list
		If list it is assumed that each elemement is array-like.
	bin_type : {'number','width','edges','equal'}, optional
		Defines how is understood the value given in bins: 'number' for the desired number of bins, 'width' for the width
		of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements (or as close
		as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if float and 'edges'
		if ndarray.
	bins : int, float, array-like or list, optional
		Gives the values for the bins, according to bin_type.
	band_stat_low : str, int, float or function, optional
		Defines how to calculate the lower limit of the error band. When passing a string it must be either one of the options
		for scipy.stats.binned_statistic(), or a string that combines 'std' with a number (e.g., '2.2std'), where to number is
		interpreted as the number of standard deviations that the limit must cover. When passing an integer or float is
		interpreted as being the percentile for the limit. When passing a function it must have the input and ouput
		characteristics required by scipy.stats.binned_statistic().
	band_stat_high : str, int, float or function, optional
		Defines how to calculate the upper limit of the error band. When passing a string it must be either one of the options
		for scipy.stats.binned_statistic(), or a string that combines 'std' with a number (e.g., '2.2std'), where to number is
		interpreted as the number of standard deviations that the limit must cover. When passing an integer or float is
		interpreted as being the percentile for the limit. When passing a function it must have the input and ouput
		characteristics required by scipy.stats.binned_statistic().
	line_stat : str, int, float or function, optional
		Defines how to calculate the centre of the error band. This centre is used to position the band when either of the band
		bounds is the standard deviation, or a multiple of that. When passing a string it must be either one of the options
		for scipy.stats.binned_statistic() When passing an integer or float is interpreted as being the percentile for the limit.
		When passing a function it must have the input and ouput characteristics required by scipy.stats.binned_statistic().
	line : boolean, optional
		If True, draw a line that follows the statistic defined in line_stat.
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
	plabel : str, optional
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
		kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, antialiasing, etc.
		A list of available `Line2D` properties can be found here: 
		https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
	
	Returns
	-------
	None
	"""
	
	import numpy as np
	from numbers import Number
	import scipy.stats as stats
	from numpy import percentile
	from functools import partial
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from splotch.base_func import axes_handler,bin_axis,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if ylog is None:
		from splotch.defaults import Params
		ylog=Params.hist1D_yaxis_log
	if bins is None:
		bins=int((len(x))**0.4)
	if 'linewidth' not in band_kw.keys():
		band_kw['linewidth']=0
	if 'alpha' not in band_kw.keys():
		band_kw['alpha']=0.4
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#band_par = {**plot_kw, **kwargs} # For Python > 3.5
	band_kw.update(kwargs)
	
	band_stat=[band_stat_low,band_stat_high]
	band_multi=np.ones(2)
	for i in range(len(band_stat)):
		if isinstance(band_stat[i],Number):
			band_stat[i]=partial(percentile,q=band_stat[i])
		elif 'std' in band_stat[i] and len(band_stat[i].replace('std',''))>0:
			band_multi[i]=float(band_stat[i].replace('std',''))
			band_stat[i]='std'
	
	if isinstance(line_stat,Number):
		line_stat=partial(percentile,q=line_stat)
	
	temp_x,bins_hist,bins_plot=bin_axis(x,bin_type,bins,log=xlog)        
	temp_y=stats.binned_statistic(temp_x,y,statistic=line_stat,bins=bins_hist)[0]
	if band_stat[0]==band_stat[1]:
		band_low,band_high=[stats.binned_statistic(temp_x,y,statistic=band_stat[0],bins=bins_hist)[0]]*2
	else:
		band_low=stats.binned_statistic(temp_x,y,statistic=band_stat[0],bins=bins_hist)[0]
		band_high=stats.binned_statistic(temp_x,y,statistic=band_stat[1],bins=bins_hist)[0]
	if band_stat[0]=='std':
		band_low=temp_y-band_multi[0]*band_low
	if band_stat[1]=='std':
		band_high=temp_y+band_multi[1]*band_high
	if ylog:
		temp_y=np.where(temp_y==0,np.nan,temp_y)
	y=np.array([temp_y[0]]+[j for j in temp_y])
	
	plt.fill_between((bins_plot[:-1]+bins_plot[1:])/2,band_low,band_high,**band_kw)
	if line:
		plt.plot((bins_plot[:-1]+bins_plot[1:])/2,temp_y,label=plabel,**line_kw)
	if plabel is not None:
		plt.legend(loc=lab_loc)
	
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

def errorbar(x,y,xerr=None,yerr=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,
				title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""Errorbar plotting function.
	
	This is a wrapper for pyplot.errorbar().
	
	Parameters
	----------
	x : array-like or list
		If list it is assumed that each elemement is array-like.
	y : array-like or list
		If list it is assumed that each elemement is array-like.
	xerr : array-like or list, optional
		Defines the length of the errobars in the x-axis. If list it is assumed that each elemement is array-like.
	yerr : array-like or list, optional
		Defines the length of the errobars in the y-axis. If list it is assumed that each elemement is array-like.
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
	plabel : str, optional
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
		kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, antialiasing, etc.
		A list of available `Line2D` properties can be found here: 
		https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
	
	Returns
	-------
	None
	"""
	
	from matplotlib.pyplot import errorbar, legend
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list:
		x=[x]
	if type(y) is not list:
		y=[y]
	if type(xerr) is not list:
		xerr=[xerr]
	if type(yerr) is not list:
		yerr=[yerr]
	L=len(x)
	if type(plabel) is not list:
		plabel=[plabel]*L

	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)

	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par = dict_splicer(plot_par,L,[1]*L)

	for i in range(L):
		errorbar(x[i],y[i],xerr=xerr[i],yerr=yerr[i],label=plabel[i],**plot_par[i])
	if any(plabel):
		legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

#Errorboxes
def errbox(x,y,xerr=None,yerr=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,boxtype='ellipse',
	title=None,xlabel=None,ylabel=None,plabel=None,grid=None,lab_loc=0,ax=None,plot_kw={},**kwargs):
	
	"""Errorbox plotting function (legacy name).
	
	This provides legacy compatibily for old code using the original name of the function.
	This function will eventually be removed, so consider switching to plots_2d.errorbox().
	"""
	
	import warnings
	
	warnings.warn('This function provides legacy support and it will be removed in the future. Use plots_2d.errorbox() instead.',FutureWarning)
	
	errorbox(x,y,xerr,yerr,xlim,ylim,xinvert,yinvert,xlog,ylog,boxtype,title,xlabel,ylabel,plabel,grid,lab_loc,ax,plot_kw,**kwargs)

def errorbox(x,y,xerr=None,yerr=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,boxtype='ellipse',
			title=None,xlabel=None,ylabel=None,plabel=None,grid=None,lab_loc=0,ax=None,plot_kw={},**kwargs):
	
	"""Errorbox plotting function.
	
	This is a wrapper around matplotlib PatchCollections with a matplotlib errorbar functionality.
	
	Parameters
	----------
	x : array-like or list
		If list it is assumed that each elemement is array-like.
	y : array-like or list
		If list it is assumed that each elemement is array-like.
	xerr : array-like or list, optional
		Defines the length of the errobars in the x-axis. If list it is assumed that each elemement is array-like.
	yerr : array-like or list, optional
		Defines the length of the errobars in the y-axis. If list it is assumed that each elemement is array-like.
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
	boxtype : str
		The type of box to plot, patch types include: ellipse | rectangle (Default: ellipse).
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	plabel : str, optional
		Sets the legend for the plot.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Patches properties.
	**kwargs: Patch properties, optional
		kwargs are used to specify matplotlib specific properties such as facecolor, linestyle, alpha, etc.
		A list of available `Patch` properties can be found here: 
		https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.Rectangle.html

	
	Returns
	-------
	None
	"""
	
	from matplotlib.pyplot import errorbar, legend
	from .base_func import axes_handler,dict_splicer,plot_finalizer

	from numpy import shape, full, array

	from matplotlib.collections import PatchCollection
	from matplotlib.patches import Ellipse, Rectangle
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list:
		x=[x]
	if type(y) is not list:
		y=[y]
	if type(xerr) is not list:
		xerr=[xerr]
	if type(yerr) is not list:
		yerr=[yerr]

	L=len(x)
	if type(plabel) is not list:
		plabel=[plabel]*L

	# Validate format of xerr and yerr
	for i in range(L):
		# x-axis errors
		if (shape(xerr[i]) == ()): # single error for all points
			xerr[i] = full((2,len(x[i])), xerr[i])
		else:
			if (len(shape(xerr[i])) == 1):
				if (shape(xerr[i])[0] == len(x[i])): # single error for each point
					xerr[i] = array([xerr[i], xerr[i]])
				elif (shape(xerr[i])[0] == 2): # separate upper and lower errors for all points
					xerr[i] = full((len(x[i]), 2), xerr[i]).T
				else:
					print('ding') # Raise exception for invalid length of points
			elif (len(shape(xerr[i])) == 2): # separate upper and lower errors for each point
				xerr[i] = array(xerr[i])
				if (shape(xerr[i])[0] != 2 or shape(xerr[i])[1] != len(x[i])):
					print('dong') # Raise exception for invalid length of points

		# y-axis errors
		if (shape(yerr[i]) == ()): # single error for all points
			yerr[i] = full((2,len(y[i])), yerr[i])
		else:
			if (len(shape(yerr[i])) == 1):
				if (shape(yerr[i])[0] == len(y[i])): # single error for each point
					yerr[i] = array([yerr[i], yerr[i]])
				elif (shape(yerr[i])[0] == 2): # separate upper and lower errors for all points
					yerr[i] = full((len(y[i]), 2), yerr[i]).T
				else:
					print('ding') # Raise exception for invalid length of points
			elif (len(shape(yerr[i])) == 2): # separate upper and lower errors for each point
				yerr[i] = array(yerr[i])
				if (shape(yerr[i])[0] != 2 or shape(yerr[i])[1] != len(y[i])):
					print('dong') # Raise exception for invalid length of points


	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)

	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	plot_par = dict_splicer(plot_par,L,[1]*L)
 
	PathColls = []
	# Loop over data points; create box/ellipse from errors at each point
	for i in range(L):
		errorboxes = []
		for xx, yy, xe, ye in zip(x[i], y[i], xerr[i].T, yerr[i].T):
			if (boxtype.lower().startswith('rect')):
				errorboxes.append( Rectangle((xx - xe[0], yy - ye[0]), xe.sum(), ye.sum()) )
			elif (boxtype.lower().startswith('ell')):
				errorboxes.append( Ellipse((xx - xe[0], yy - ye[0]), xe.sum(), ye.sum()) )
			else:
				print('dang')

		# Create and add patch collection with specified colour/alpha
		pc = PatchCollection(errorboxes, **plot_par[i])
		ax.add_collection(pc)

	# for i in range(L):
	# 	errorbar(x[i],y[i],xerr=xerr[i],yerr=yerr[i],label=plabel[i],**plot_par[i])
	if any(plabel):
		legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)


# Histogram and 2D binned statistics
def hist2D(x,y,bin_type=None,bins=None,dens=True,scale=None,c=None,cstat=None,xlim=None,ylim=None,clim=[None,None],nmin=0, 
			xinvert=False,yinvert=False,cbar_invert=False,xlog=False,ylog=False,clog=None,title=None,xlabel=None,
			ylabel=None,clabel=None,lab_loc=0,ax=None,grid=None,output=None,plot_kw={},**kwargs):
	
	"""2D histogram function.
	
	Parameters
	----------
	x : array-like
		Position of data points in the x axis.
	y : array-like
		Position of data points in the y axis.
	bin_type : {'number','width','edges','equal'}, optional
		Defines how is understood the value given in bins: 'number' for givinf the desired number of bins, 'width' for
		the width of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements
		(or as close as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if
		float and 'edges' if ndarray.
	bins : int, float, array-like or list, optional
		Gives the values for the bins, according to bin_type.
	dens : bool or list, optional
		If false the histogram returns raw counts.
	scale : float or list, optional
		Scaling of the data counts.
	c : array-like, optional
		If a valid argument is given in cstat, defines the value used for the binned statistics.
	cstat : str or function, optional
		Must be one of the valid str arguments for the statistics variable in scipy.stats.binned_statistic_2d
		('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or a function that takes a 1D array and outputs an integer
		 or float.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	clim : list, optional
		Defines the limits of the colour map ranges, it must contain two elements (lower and higer limits).
	nmin : int, optional (default: 0)
		The minimum number of points required in a bin in order to be plotted.
	xinvert : bool, optional
		If true inverts the x-axis.
	yinvert : bool, optional
		If true inverts the y-axis.
	cbar_invert : bool, optional
		If True inverts the direction of the colour bar (not the colour map).
	xlog : bool, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool, optional
		If True the scale of the x-axis is logarithmic.
	clog : bool, optional
		If True, the colour map is changed from linear to logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	clabel : str, optional
		Sets the legend for the colour axis.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	output : boolean, optional
		If True, returns the edges and values of the histogram.
	plot_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib scatter function.
		Parameters will be overwritten if also given implicitly as a **kwarg.
	**kwargs : pcolormesh properties, optional
		kwargs are used to specify matplotlib specific properties such as cmap, norm, edgecolors etc.
		https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html

	Returns
	-------
	n : array
		The values of the histogram. Only provided if output is True.
	bin_edges_x : array
		The bin edges for the x axis. Only provided if output is True.
	bin_edges_y : array
		The bin edges for the y axis. Only provided if output is True.
	"""
	
	from numpy import nan
	from matplotlib.colors import LogNorm
	from matplotlib.pyplot import pcolormesh, colorbar
	from .base_func import axes_handler,basehist2D,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(bin_type) is not list:
		bin_type=[bin_type]*2
	if type(bins) not in [list,tuple]:
		if bins is None:
			bins=int((len(x))**0.4)
		bins=[bins]*2
	X,Y,Z=basehist2D(x,y,c,bin_type,bins,scale,dens,cstat,xlog,ylog)
	_,_,counts = basehist2D(x,y,c,bin_type,bins,scale,dens,'count',xlog,ylog) # Also get counts for number threshold cut

	# Cut bins which do not meet the number count threshold
	Z[counts<nmin] = nan
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)
	
	if None in (clog,output):
		from .defaults import Params
		if clog is None:
			clog=Params.hist2D_caxis_log
		if output is None:
			output=Params.hist2D_output
	if clog:
		pcolormesh(X,Y,Z.T,norm=LogNorm(vmin=clim[0],vmax=clim[1],clip=True),**plot_par)
	else:
		if cstat is None:
			Z[Z==0]=nan
		pcolormesh(X,Y,Z.T,vmin=clim[0],vmax=clim[1],**plot_par)
	if clabel is not None:
		cbar=colorbar()
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)
	if output:
		return(Z.T,X,Y)

# Image
def img(im,x=None,y=None,xlim=None,ylim=None,clim=[None,None],cmin=0,xinvert=False,yinvert=False,cbar_invert=False,clog=None,
		title=None,xlabel=None,ylabel=None,clabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""2D pixel-based image plotting function.
	
	Parameters
	----------
	im : array-like
		Value for each pixel in an x-y 2D array, where the first dimension is the x-position and the second is
		the y-position.
	x : array-like, optional
		Position of data points in the x axis.
	y : array-like, optional
		Position of data points in the y axis.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	clim : list, optional
		Defines the limits of the colour map ranges, it must contain two elements (lower and higer limits).
	clog : bool, optional
		If True, the colour map is changed from linear to logarithmic.
	xinvert : bool, optional
		If true inverts the x-axis.
	yinvert : bool, optional
		If true inverts the y-axis.
	cbar_invert : bool, optional
		If True inverts the direction of the colour bar (not the colour map).
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	clabel : str, optional
		Sets the legend for the colour axis.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	plot_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib pcolormesh function.
		Parameters will be overwritten if also given implicitly as a **kwarg.
	**kwargs : pcolormesh properties, optional
		kwargs are used to specify matplotlib specific properties such as `cmap`, `marker`, `norm`, etc.
		A list of available `pcolormesh` properties can be found here:
		https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html
	
	Returns
	-------
	None
	"""
	
	from numpy import arange, meshgrid
	from matplotlib.colors import LogNorm
	from matplotlib.pyplot import pcolormesh, colorbar
	from .base_func import axes_handler,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if x is None:
		x=arange(len(im[:,0])+1)
	if y is None:
		y=arange(len(im[0,:])+1)
	if clog is None:
		from .defaults import Params
		clog=Params.img_caxis_log

	X, Y = meshgrid(x, y)

	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs) 

	if clog:
		pcolormesh(X,Y,im.T,norm=LogNorm(vmin=clim[0],vmax=clim[1],clip=True),**plot_par)
	else:
		pcolormesh(X,Y,im.T,vmin=clim[0],vmax=clim[1],**plot_par)
	if clabel is not None:
		cbar=colorbar()
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()
	plot_finalizer(False,False,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

# Scatter
def scatter(x,y,c=None,xlim=None,ylim=None,clim=None,xinvert=False,yinvert=False,cbar_invert=False,xlog=False,ylog=False,title=None,
			xlabel=None,ylabel=None,clabel=None,plabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""2D pixel-based image plotting function.
	
	Parameters
	----------
	x : array-like or list
		Position of data points in the x-axis.
	y : array-like or list
		Position of data points in the y-axis.
	c : array-like or list, optional
		Value of data points in the z-axis (colour-axis).
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	clim : tuple-like, optional
		Defines the limits of the colour-axis, it must contain two elements (lower and higer limits).
		Functions equivalently to the `vmin, vmax` arguments used by `colors.Normalize`. If both are given,
		clim` takes priority.
	xinvert : bool, optional
		If true inverts the x-axis.
	yinvert : bool, optional
		If true inverts the y-axis.
	cbar_invert : bool, optional
		If True inverts the direction of the colour bar (not the colour map).
	xlog : bool, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool, optional
		If True the scale of the x-axis is logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	clabel : str, optional
		Sets the legend for the colour axis.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	plot_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib scatter function.
		Parameters will be overwritten if also given implicitly as a **kwarg.
	**kwargs : Collection properties, optional
		kwargs are used to specify matplotlib specific properties such as cmap, marker, norm, etc.
		A list of available `Collection` properties can be found here:
		https://matplotlib.org/3.1.0/api/collections_api.html#matplotlib.collections.Collection
	
	Returns
	-------
	paths
		A list of PathCollection objects representing the plotted data.
	"""

	from numpy import shape
	from matplotlib.pyplot import scatter, colorbar, legend
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list or len(shape(x))==1:
		x=[x]
	if type(y) is not list or len(shape(y))==1:
		y=[y]
	if type(c) is not list or len(shape(c))==1:
		c=[c]
	L=len(x)
	if type(plabel) is not list:
		plabel=[plabel]*L

	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)

	# Insert clim as vmin, vmax into **kwargs dictionary, if given.
	if (clim != None):
		try:
			_ = (e for e in clim)
			if (len(clim) == 2):
				plot_par['vmin'] = clim[0]
				plot_par['vmax'] = clim[1]
			else:
				raise TypeError("`clim` must be of iterable type and have two values only.")
		except (TypeError):
			raise TypeError("`clim` must be of iterable type and have two values only.")

	# Create 'L' number of plot kwarg dictionaries to parse into each scatter call
	plot_par=dict_splicer(plot_par,L,[len(i) for i in x])

	paths = []
	for i in range(L):
		p = scatter(x[i],y[i],c=c[i],label=plabel[i],**plot_par[i])
		paths.append(p)
	if clabel is not None:
		cbar = colorbar()
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()
	if any(plabel):
		legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

	return paths[0] if len(paths) == 1 else paths

def sector(r,theta,rlim=(0.0,1.0),thetalim=(0.0,360.0),clim=None,rotate=0.0,rlabel="",thetalabel="",clabel=None,rstep=None,
			thetastep=15.0,rticks='auto',thetaticks='auto',cbar_invert=False,fig=None,plot_kw={},**kwargs):
	
	""" Sector Plot function
	
	Plots a sector plot (a.k.a "pizza plot") based on data with one radial axis and an angular axis
	
	Parameters
	----------
	r : array-like or list
		Radial axis data.
	theta : array-like or list
		Angular axis data (degrees).
	rlim : tuple-like, optional
		The lower and upper limits for the radial axis (degrees).
	thetalim : tuple-like, optional
		The lower and upper limits for the angular axis (degrees).
	clim : tuple-like, optional
		Defines the limits of the colour-axis, it must contain two elements (lower and higer limits).
		Functions equivalently to the `vmin, vmax` arguments used by `colors.Normalize`. If both are given,
		clim` takes priority.
	rotate : float, optional
		By how many degrees (clockwise) to rotate the entire plot (valid values in [-180, 180]).
	rlabel : str, optional
		Sets the label of the r-axis.
	thetalabel : str, optional
		Sets the label of the theta-axis.
	clabel : str, optional
		Sets the legend for the colour axis.
	rstep : float, optional
		Sets the step size of r ticks.
	thetastep : float, optional, default: 15.0
		Sets the step size of theta ticks (degrees).
	rticks : 'auto', or ticker
		* Not implement *
	thetaticks : 'auto', or ticker
		* Not implement *
	cbar_invert : bool, optional
		If True inverts the direction of the colour bar (not the colour map).
	fig : pyplot.Figure, optional
		Use the given figure to make the plot, defaults to the current figure.
	plot_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib scatter function.
		Parameters will be overwritten if also given implicitly in **kwargs.
	**kwargs : Collection properties, optional
		kwargs are used to specify matplotlib specific properties such as cmap, marker, norm, etc.
		A list of available `Collection` properties can be found here:
		https://matplotlib.org/3.1.0/api/collections_api.html#matplotlib.collections.Collection
	
	Returns
	-------
	ax : The pyplot.Axes object created for the sector plot.
	"""
	
	from matplotlib.transforms import Affine2D
	from matplotlib.projections.polar import PolarAxes
	from matplotlib.pyplot import gcf, colorbar, legend 
	
	from mpl_toolkits.axisartist import floating_axes
	from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator, DictFormatter)
	import mpl_toolkits.axisartist.angle_helper as angle_helper
	
	from numpy import array, linspace, arange, shape, sqrt, floor, round, degrees, radians, pi
	
	if (fig == None):
		fig = gcf()
	
	# rotate a bit for better orientation
	trans_rotate = Affine2D().translate(0.0, 0)
	
	# scale degree to radians
	trans_scale = Affine2D().scale(pi/180.0, 1.)
	trans = trans_rotate + trans_scale + PolarAxes.PolarTransform()
	
	# Get theta ticks
	#if (thetaticks == 'auto'):
	thetaticks = arange(*radians(array(thetalim)-rotate),step=radians(thetastep))
	theta_gridloc = FixedLocator(thetaticks[thetaticks/(2*pi) < 1])
	theta_tickfmtr = DictFormatter(dict(zip(thetaticks,[f"{(round(degrees(tck)+rotate)):g}" for tck in thetaticks])))
	
	#tick_fmtr = DictFormatter(dict(angle_ticks))
	#tick_fmtr = angle_helper.Formatter()
	
	if (rstep == None):
		rstep = 0.5
	
	r_gridloc = FixedLocator(arange(rlim[0],rlim[1],step=rstep))
	
	grid = floating_axes.GridHelperCurveLinear(
		PolarAxes.PolarTransform(),
		extremes=(*radians(array(thetalim)-rotate), *rlim),
		grid_locator1=theta_gridloc,
		grid_locator2=r_gridloc,
		tick_formatter1=theta_tickfmtr,
		tick_formatter2=None,
	)
	
	ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid)
	fig.add_subplot(ax)
	
	# tick references
	thetadir_ref = ['top','right','bottom','left']
	rdir_ref = ['bottom','left','top','right']
	
	# adjust axes directions
	ax.axis["left"].set_axis_direction('bottom') # Radius axis (displayed)
	ax.axis["right"].set_axis_direction('top') # Radius axis (hidden)
	ax.axis["top"].set_axis_direction('bottom') # Theta axis (outer)
	ax.axis["bottom"].set_axis_direction('top') # Theta axis (inner)
	
	# Top theta axis
	ax.axis["top"].toggle(ticklabels=True, label=True)
	ax.axis["top"].major_ticklabels.set_axis_direction(thetadir_ref[(int(rotate)//90)%4])
	ax.axis["top"].label.set_axis_direction(thetadir_ref[(int(rotate)//90)%4])
	
	# Bottom theta axis
	ax.axis["bottom"].set_visible(False if rlim[0] < (rlim[1]-rlim[0])/3 else True)
	ax.axis["bottom"].major_ticklabels.set_axis_direction(thetadir_ref[(int(rotate)//90+2)%4])
    
	# Visible radius axis    
	ax.axis["left"].major_ticklabels.set_axis_direction(rdir_ref[(int(rotate)//90)%4])
	ax.axis["left"].label.set_axis_direction(rdir_ref[(int(rotate)//90)%4])
    
	# Labels
	ax.axis["left"].label.set_text(rlabel)
	ax.axis["top"].label.set_text(thetalabel)
	
	# create a parasite axes whose transData in RA, cz
	sector_ax = ax.get_aux_axes(trans)
	
	# This has a side effect that the patch is drawn twice, and possibly over some other
	# artists. So, we decrease the zorder a bit to prevent this. 
	sector_ax.patch = ax.patch  
	sector_ax.patch.zorder = 0.9
	
	
	L = shape(theta)[0] if len(shape(theta)) > 1 else 1
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)

	# Insert clim as vmin, vmax into **kwargs dictionary, if given.
	if (clim != None):
		try:
			_ = (e for e in clim)
			if (len(clim) == 2):
				plot_par['vmin'] = clim[0]
				plot_par['vmax'] = clim[1]
			else:
				raise TypeError("`clim` must be of iterable type and have two values only.")
		except (TypeError):
			raise TypeError("`clim` must be of iterable type and have two values only.")
	
	# Create 'L' number of plot kwarg dictionaries to parse into each plot call
	#plot_par = dict_splicer(plot_par,L,[1]*L)
	
	if (L == 1):
		sctr = sector_ax.scatter(theta-rotate, r, **plot_par)
	else:
		for ii in range(L):
			sctr = sector_ax.scatter(theta[ii]-rotate, r[ii],**plot_par[ii])
	
	if clabel is not None:
		cbar = colorbar(sctr)
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()

	return sector_ax

# Contours encircling the densest region down to a certain percentage 
def sigma_cont(x,y,percent=None,bin_type=None,bins=None,c=None,cmap=None,xlim=None,ylim=None,clim=None,xinvert=False,yinvert=False,
				cbar_invert=False,s=None,xlog=False,ylog=False,title=None,xlabel=None,ylabel=None,clabel=None,lab_loc=0,ax=None,
				grid=None,output=None,plot_kw={},**kwargs):
	
	"""Contour function, encircling the highest density regions that contain the given percentages of the sample (legacy name).
	
	This provides legacy compatibily for old code using the original name of the function.
	This function will eventually be removed, so consider switching to plots_2d.contourp().
	"""
	
	import warnings
	
	warnings.warn('This function provides legacy support and it will be removed in the future. Use plots_2d.contourp() instead.',FutureWarning)
	
	contourp(x,y,percent,bin_type,bins,c,cmap,xlim,ylim,clim,xinvert,yinvert,cbar_invert,s,xlog,ylog,title,xlabel,ylabel,
				clabel,lab_loc,ax,grid,output,plot_kw,**kwargs)

def contourp(x,y,percent=None,bin_type=None,bins=None,c=None,cmap=None,xlim=None,ylim=None,clim=None,xinvert=False,yinvert=False,
				cbar_invert=False,s=None,xlog=False,ylog=False,title=None,xlabel=None,ylabel=None,clabel=None,lab_loc=0,ax=None,
				grid=None,output=None,plot_kw={},**kwargs):
	
	"""Contour function, encircling the highest density regions that contain the given percentages of the sample.
	
	Parameters
	----------
	x : array-like
		Position of data points in the x axis.
	y : array-like
		Position of data points in the y axis.
	percent : float or list, optional.
		The percentages of the sample that the contours encircle.
	bin_type : {'number','width','edges','equal'}, optional
		Defines how is understood the value given in bins: 'number' for givinf the desired number of bins, 'width' for
		the width of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements
		(or as close as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if
		float and 'edges' if ndarray.
	bins : int, float, array-like or list, optional
		Gives the values for the bins, according to bin_type.
	c : str, float or list, optional
		The colours of the contours. If float or list of floats, they must be in the range [0,1], as the colours are
		taken from the given colour map.
	cmap : str, optional
		The colour map to be used, viridis by default.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
	clim : list, optional
		Defines the limits of the colour map ranges, it must contain two elements (lower and higer limits).
	xinvert : bool, optional
		If true inverts the x-axis.
	yinvert : bool, optional
		If true inverts the y-axis.
	cbar_invert : bool, optional
		If True inverts the direction of the colour bar (not the colour map).
	s= str or list, optional.
		Defines the linestyle of the contours.
	xlog : bool, optional
		If True the scale of the x-axis is logarithmic.
	ylog : bool, optional
		If True the scale of the x-axis is logarithmic.
	title : str, optional
		Sets the title of the plot
	xlabel : str, optional
		Sets the label of the x-axis.
	ylabel : str, optional
		Sets the label of the y-axis.
	clabel : str, optional
		Sets the legend for the colour axis.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	grid : boolean, optional
		If not given defaults to the value defined in splotch.Params.
	output : boolean, optional
		If True, returns the edges and values of the underlying histogram plus the levels of the contours.
	plot_kw : dict, optional
		Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are QuadContourSet properties.
	**kwargs: QuadContourSet properties, optional
		kwargs are used to specify matplotlib specific properties such as cmap, linewidths, hatches, etc.
		The list of available properties can be found here: 
		https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.contour.html
	
	Returns
	-------
	bin_edges_x : array
		The bin edges for the x axis.
	bin_edges_y : array
		The bin edges for the y axis.
	n : array
		The values of the underlying histogram.
	l : array
		The levels for the contours.
	"""
	
	from matplotlib.cm import get_cmap
	from numpy import array,linspace, round
	from matplotlib.pyplot import gca, contour, legend
	from .base_func import axes_handler,basehist2D,percent_finder,plot_finalizer,dict_splicer
	
	if None in (percent,cmap,clim,s,output):
		from .defaults import Params
		if percent is None:
			percent=Params.sigcont_percent
		if cmap is None:
			cmap=Params.sigcont_cmap
		if clim is None:
			clim=Params.sigcont_clim
		if s is None:
			s=Params.sigcont_linestyle
		if output is None:
			output=Params.sigcont_output
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(percent) is not list:
		percent=[percent]
	if type(bin_type) is not list:
		bin_type=[bin_type]*2
	if type(bins) is not list:
		if bins is None:
			bins=int((len(x))**0.4)
		bins=[bins]*2
	if output is None:
		from .defaults import Params
		output=Params.hist2D_output
	
	cmap=get_cmap(cmap)
	X,Y,Z=basehist2D(x,y,c,bin_type,bins,None,None,None,xlog,ylog)
	X=(X[:-1]+X[1:])/2
	Y=(Y[:-1]+Y[1:])/2
	CS=[]
	if c is None:
		if len(percent)<4:
			col_ax=gca()
			l=col_ax.plot([1,2,3])
			c=[l[0].get_color()]*len(percent)
			l.pop(0).remove()
		else:
			if len(s)<4:
				s=['solid']*len(percent)
			c=cmap(linspace(clim[0],clim[1],len(percent)))
	else:
		if type(c) is str:
			c=[c]*len(percent)
		else:
			c=cmap(c)
			s=['solid']*len(percent)
	if type(clabel) is not list:
		if clabel is None:
			clabel=[str(round(p,1))+'%' for p in percent]
		else:
			clabel= [clabel] + [None]*(len(percent)-1)
	if type(s) is not list:
		s=[s]*len(percent)
	
	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)
	
	# Create 'L' number of plot kwarg dictionaries to parse into each scatter call
	plot_par=dict_splicer(plot_par,len(percent),[1]*len(percent))
	
	level=[]
	for i in range(len(percent)):
		level.append(percent_finder(Z,percent[i]/100))
		CS.append(contour(X,Y,Z.T,levels=[level[i]],colors=[c[i],],linestyles=s[i],**plot_par[i]))
		if clabel[0] is not None:
			CS[i].collections[0].set_label(clabel[i])
	if clabel[0] is not None:
		legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)
	if output:
		return(X,Y,Z.T,array(level))

