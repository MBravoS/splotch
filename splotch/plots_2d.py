#### Definition of all wrappers for 2D plotting

#Errorbars
def errbar(x,y,xerr=None,yerr=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,title=None,
			xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	import numpy as np
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
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
		plt.errorbar(x[i],y[i],xerr=xerr[i],yerr=yerr[i],label=plabel[i],**plot_par[i])
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

# Histogram and 2D binned statistics
def hist2D(x,y,bin_type=None,bins=None,dens=True,scale=None,c=None,cstat=None,xlim=None,ylim=None,clim=[None,None],
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
	
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,base_hist2D,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(bin_type) is not list:
		bin_type=[bin_type]*2
	if type(bins) is not list:
		if bins is None:
			bins=int((len(x))**0.4)
		bins=[bins]*2
	X,Y,Z=base_hist2D(x,y,c,bin_type,bins,scale,dens,cstat,xlog,ylog)

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
		plt.pcolormesh(X,Y,Z.T,norm=clr.LogNorm(vmin=clim[0],vmax=clim[1],clip=True),**plot_par)
	else:
		if cstat is None:
			Z[Z==0]=np.nan
		plt.pcolormesh(X,Y,Z.T,vmin=clim[0],vmax=clim[1],**plot_par)
	if clabel is not None:
		cbar=plt.colorbar()
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)
	if output:
		return(Z.T,X,Y)

# Image
def img(im,x=None,y=None,xlim=None,ylim=None,clim=[None,None],xinvert=False,yinvert=False,cbar_invert=False,clog=None,
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
	
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if x is None:
		x=np.arange(len(im[:,0])+1)
	if y is None:
		y=np.arange(len(im[0,:])+1)
	if clog is None:
		from .defaults import Params
		clog=Params.img_caxis_log

	X, Y = np.meshgrid(x, y)

	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs) 

	if clog:
		plt.pcolormesh(X,Y,im.T,norm=clr.LogNorm(vmin=clim[0],vmax=clim[1],clip=True),**plot_par)
	else:
		plt.pcolormesh(X,Y,im.T,vmin=clim[0],vmax=clim[1],**plot_par)
	if clabel is not None:
		cbar=plt.colorbar()
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()
	plot_finalizer(False,False,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

# Scatter
def scatter(x,y,c=None,xlim=None,ylim=None,xinvert=False,yinvert=False,cbar_invert=False,xlog=False,ylog=False,title=None,
			xlabel=None,ylabel=None,clabel=None,plabel=None,lab_loc=0,ax=None,grid=None,plot_kw={},**kwargs):
	
	"""2D pixel-based image plotting function.
	
	Parameters
	----------
	x : array-like or list
		Position of data points in the x axis.
	y : array-like or list
		Position of data points in the y axis.
	c : array-like or list, optional
		Value of data points in the z axis (colour axis).
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
		Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
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
	None
	"""
	
	import numpy as np
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list or len(np.shape(x))==1:
		x=[x]
	if type(y) is not list or len(np.shape(y))==1:
		y=[y]
	if type(c) is not list or len(np.shape(c))==1:
		c=[c]
	L=len(x)
	if type(plabel) is not list:
		plabel=[plabel]*L

	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)

	# Create 'L' number of plot kwarg dictionaries to parse into each scatter call
	plot_par=dict_splicer(plot_par,L,[len(i) for i in x])

	for i in range(L):
		plt.scatter(x[i],y[i],c=c[i],label=plabel[i],**plot_par[i])
	if clabel is not None:
		cbar=plt.colorbar()
		cbar.set_label(clabel)
		if cbar_invert:
			cbar.ax.invert_yaxis()
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)

# Contours encircling the densest part down to a certain percetange 
def sigma_cont(x,y,percent=[68.27,95.45],bin_type=None,bins=None,c=None,cmap='viridis',xlim=None,ylim=None,
				clim=[0.33,0.67],xinvert=False,yinvert=False,cbar_invert=False,s=['solid','dashed','dotted'],xlog=False,
				ylog=False,title=None,xlabel=None,ylabel=None,clabel=None,lab_loc=0,ax=None,grid=None,output=None,plot_kw={},**kwargs):
	
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
		If True, returns the edges and values of the histogram.
	plot_kw : dict, optional
		Explicit dictionary of kwargs to be parsed to matplotlib scatter function.
		Parameters will be overwritten if also given implicitly as a **kwarg.
	
	Returns
	-------
	bin_edges_x : array
		The bin edges for the x axis.
	bin_edges_y : array
		The bin edges for the y axis.
	n : array
		The values of the underlying histogram.
	"""
	
	import numpy as np
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,base_hist2D,percent_finder,plot_finalizer
	
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
	cmap=cm.get_cmap(cmap)
	X,Y,Z=base_hist2D(x,y,c,bin_type,bins,None,None,None,xlog,ylog)
	X=(X[:-1]+X[1:])/2
	Y=(Y[:-1]+Y[1:])/2
	CS=[]
	if c is None:
		if len(percent)<4:
			col_ax=plt.gca()
			l=col_ax.plot([1,2,3])
			c=[l[0].get_color()]*len(percent)
			l.pop(0).remove()
		else:
			if len(s)<4:
				s=['solid']*len(percent)
			c=cmap(np.linspace(clim[0],clim[1],len(percent)))
	else:
		if type(c) is str:
			c=[c]*len(percent)
		else:
			c=cmap(c)
			s=['solid']*len(p)
	if type(clabel) is not list:
		if clabel is None:
			clabel=[str(np.round(p,1))+'%' for p in percent]
		else:
			clabel=[clabel]*len(x)

	# Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
	#plot_par = {**plot_kw, **kwargs} # For Python > 3.5
	plot_par = plot_kw.copy()
	plot_par.update(kwargs)


	for i in range(len(percent)):
		level=[percent_finder(Z,percent[i]/100)]
		CS.append(plt.contour(X,Y,Z.T,levels=level,colors=[c[i],],linewidths=1.5,linestyles=s[i]))
		if clabel[0] is not None:
			CS[i].collections[0].set_label(clabel[i])
	if clabel[0] is not None:
		plt.legend(loc=lab_loc)
	plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid)
	if ax is not None:
		old_axes=axes_handler(old_axes)
	if output:
		return(Z.T,X,Y)
