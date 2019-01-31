#### Definition of all wrappers for 1D plotting

#Histogram
def hist(data,bin_num=None,dens=True,norm=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=True,
			title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,plot_par={},multi=False):
	
	"""1D histogram function.
	
	The plotting is done with pyplot.plot(), so histograms are shown with interpolated curves instead of the
	more common stepwise curve. For this reason splotch.histstep is a better choice for small datasets. 
	
	Parameters
	----------
	data : array-like or list
		If list it is assumed that each elemement is array-like.
	bin_num : int or list, optional
		Number of bins.
	dens :  bool or list, optional
		If false the histogram returns raw counts.
	norm : float or list, optional
		Normalization of the counts.
	xlim : tuple-like, optional
		Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
	ylim : tuple-like, optional
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
	plabel : str, optional
		Sets the legend for the plot.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	plot_par : dict, optional
		Passes the given dictionary as a kwarg to the plotting function.
	multi : bool, optional
		If True, holds the application of x/ylog, x/yinvert and grid, to avoid duplication.
	
	Returns
	-------
	None
	"""
	
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(data) is not list:
		data=[data]
	L=len(data)
	if bin_num is None:
		bin_num=[int((len(d))**0.4) for d in data]
	if type(bin_num) is not list:
		bin_num=[bin_num+1]*L
	if type(dens) is not list:
		dens=[dens]*L
	if type(norm) is not list:
		norm=[norm]*L
	if type(plabel) is not list:
		plabel=[plabel]*L
	plot_par=dict_splicer(plot_par,L)
	for i in range(L):
		if xlog:
			bins=np.logspace(np.log10(np.nanmin(data[i])),np.log10(np.nanmax(data[i])),num=bin_num[i])
		else:
			bins=np.linspace(np.nanmin(data[i]),np.nanmax(data[i]),num=bin_num[i])
		y,x=np.histogram(data[i],bins=bins,density=dens[i])
		if dens[i]:
			if norm[i]:
				y*=1.0*len(data[i])/norm[i]
		plt.plot((x[0:-1]+x[1:])/2,y,label=plabel[i],**plot_par[i])
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(old_axes)

#Step histogram
def histstep(data,bin_num=None,dens=True,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=True,
			title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,plot_par={},multi=False):
	
	"""'Clasic' 1D histogram function.
	
	This function is designed around pyplot.hist(), so it lacks the functionality to use arbitraty y-axis
	normalisation of splotch.hist().
	It is better choice for small datasets, as it plots with stepwise curves, instead of the interpolated
	ones of splotch.hist().
	
	Parameters
	----------
	data : array-like or list
		If list it is assumed that each elemement is array-like.
	bin_num : int or list, optional
		Number of bins.
	dens :  bool or list, optional
		If false the histogram returns raw counts.
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
	plot_par : dict, optional
		Passes the given dictionary as a kwarg to the plotting function.
	multi : bool, optional
		If True, holds the application of x/ylog, x/yinvert and grid, to avoid duplication.
	
	Returns
	-------
	None
	"""
	
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(data) is not list:
		data=[data]
	L=len(data)
	if bin_num is None:
		bin_num=[int((len(d))**0.4) for d in data]
	if type(bin_num) is not list:
		bin_num=[bin_num+1]*L
	if type(plabel) is not list:
		plabel=[plabel]*L
	plot_par=dict_splicer(plot_par,L)
	for i in range(L):
		if xlog:
			bins=np.logspace(np.log10(np.nanmin(data[i])),np.log10(np.nanmax(data[i])),num=bin_num[i])
		else:
			if np.nanmin(data[i])==np.nanmax(data[i]):
				bins=np.linspace(np.nanmin(data[i])-0.5,np.nanmax(data[i])+0.5,num=bin_num[i])
			else:
				bins=np.linspace(np.nanmin(data[i]),np.nanmax(data[i]),num=bin_num[i])
		plt.hist(data[i],bins=bins,density=dens,label=plabel[i],**plot_par[i])
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(old_axes)

# Generalized lines   NEEDS TO BE REDESIGNED
#def line(x,y,n=10,xlim=None,ylim=None,xinvert=False,yinvert=False,cinvert=False,xlog=False,ylog=False,title=None,
#			xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,plot_par={},multi=False):
#	import numpy as np
#	import matplotlib.pyplot as plt
#	from .base_func import axes_handler,plot_finalizer
#	
#	if ax is not None:
#		old_axes=axes_handler(ax)
#	if xlog:
#		x=np.logspace(np.log10(x[0]),np.log10(x[1]),num=n)
#	else:
#		x=np.linspace(x[0],x[1],num=n)
#	if ylog:
#		y=np.logspace(np.log10(y[0]),np.log10(y[1]),num=n)
#	else:
#		y=np.linspace(y[0],y[1],num=n)
#	plt.plot(x,y,color=c,linestyle=line_style,label=plabel,**plot_par)
#	if plabel is not None:
#		plt.legend(loc=lab_loc)
#	if not multi:
#		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
#	if ax is not None:
#		old_axes=axes_handler(ax)

#Plots
def plot(x,y=None,xlim=None,ylim=None,xinvert=False,yinvert=False,xlog=False,ylog=False,title=None,xlabel=None,ylabel=None,
			plabel=None,lab_loc=0,ax=None,plot_par={},multi=False):
	import numpy as np
	import matplotlib.pyplot as plt
	from .base_func import axes_handler,dict_splicer,plot_finalizer
	
	"""Base plotting function.
	
	This is a wrapper for pyplot.plot().
	
	Parameters
	----------
	x : array-like or list
		If list it is assumed that each elemement is array-like. If y is not given, the given values pass to y and a
		numpy array is generated with numpy.arange() for the x values.
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
	plabel : str, optional
		Sets the legend for the plot.
	lab_loc : int, optional
		Defines the position of the legend
	ax : pyplot.Axes, optional
		Use the given axes to make the plot, defaults to the current axes.
	plot_par : dict, optional
		Passes the given dictionary as a kwarg to the plotting function.
	multi : bool, optional
		If True, holds the application of x/ylog, x/yinvert and grid, to avoid duplication.
	
	Returns
	-------
	None
	"""
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list:
		x=[x]
	L=len(x)
	if y is None:
		y=x
		x=[np.arange(len(x[i])) for i in range(L)]
	else:
		if type(y) is not list:
			y=[y]
	if type(plabel) is not list:
		plabel=[plabel]*L
	plot_par=dict_splicer(plot_par,L)
	for i in range(L):
		plt.plot(x[i],y[i],label=plabel[i],**plot_par[i])
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(old_axes)

