def axes_handler(new_axis):
	"""New axis handler
	
	This function is a base-level function used by most other plotting functions to set the current Axes instance
	to ax and returns the old Axes instance to be later reverted.
	
	Parameters
	----------
	new_axis : Axes object
		The new Axes instance to be set as the current axis.
	
	Returns
	-------
	curr_axis : Axes object
		The previous Axes instance.
	"""
	from matplotlib.pyplot import gca,sca
	
	curr_axis=gca()
	sca(new_axis)
	return(curr_axis)

def base_hist2D(x,y,c,bin_type,bin_num,norm,dens,cstat,xlog,ylog):
	"""2D histogram base calculation
	
	This function is a base-level function used by plots_2d.hist2D() and plots_2d.sigma_cont() to calculate the
	underlying histogram.
	
	Parameters
	----------
	x : ndarray
		Position of data points in the x axis.
	y : ndarray
		Position of data points in the y axis.
	c : ndarray
		If a valid argument is given in cstat, defines the value used for the binned statistics.
	bin_type : {'number','width','edges','equal'}
		Defines how is understood the value given in bins: 'number' for givinf the desired number of bins, 'width' for
		the width of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements
		(or as close as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if
		float and 'edges' if ndarray.
	bin_num : int, float or array-like
		Gives the values for the bins, according to bin_type.
	norm : float
		Normalization of the counts.
	dens : bool
		If false the histogram returns raw counts.
	cstat : str or function
		Must be one of the valid str arguments for the statistics variable in scipy.stats.binned_statistic_2d
		('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or a function that takes a 1D array and outputs an integer
		 or float.
	xlog : bool
		If True the scale of the x-axis is logarithmic.
	ylog : bool
		If True the scale of the x-axis is logarithmic.
	
	Returns
	-------
	x_bins_plot : ndarray
		The bin edges on the x-axis.
	y_bins_plot : ndarray
		The bin edges on the y-axis.
	Z : ndarray
		The value of each bin.
	"""
	from numpy import histogram2d
	from scipy.stats import binned_statistic_2d
	
	x_temp,x_bins_hist,x_bins_plot=binned_axis(x,bin_type[0],bin_num[0],log=xlog)
	y_temp,y_bins_hist,y_bins_plot=binned_axis(y,bin_type[1],bin_num[1],log=ylog)
	if cstat:
		Z=binned_statistic_2d(x_temp,y_temp,c,statistic=cstat,bins=[x_bins_hist,y_bins_hist])[0]
	else:
		Z=histogram2d(x_temp,y_temp,bins=[x_bins_hist,y_bins_hist],density=dens)[0]
		if dens and norm:
			Z*=1.0*len(x)/norm
	return(x_bins_plot,y_bins_plot,Z)

def binned_axis(data,btype,bins,log=False):
	"""Bin construction for histograms
	
	This function is a base-level function used by all histogram-related functions to construct the bins.
	
	Parameters
	----------
	data : ndarray
		Data to be binned.
	bin_type : {'number','width','edges','equal'}
		Defines how is understood the value given in bins: 'number' for givinf the desired number of bins, 'width' for
		the width of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements
		(or as close as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if
		float and 'edges' if ndarray.
	bins : int, float, array-like or list
		Gives the values for the bins, according to bin_type.
	log: bool, optional
		If True, the bins are constructed in logarithmic space and the logarithm of the data is returned. 
	
	Returns
	-------
	data : ndarray
		The data for the histogram.
	hist_bins : ndarray
		The bin edges.
	plot_bins: ndarray
		The bin centers.
	
	"""
	from numpy import linspace,nanmax,nanmin
	
	def N(d,b):
		if nanmin(d)==nanmax(d):
			h=linspace(nanmin(d)-0.5,nanmax(d)+0.5,num=b)
		else:
			h=linspace(nanmin(d),nanmax(d),num=b+1)
		return(h)
	
	def W(d,b):
		from math import ceil
		if nanmin(d)==nanmax(d):
			h=array([nanmin(d)-b/2,nanmax(d)+b/2])
		else:
			L=ceil((nanmax(d)-nanmin(d))/b)
			h=nanmin(d)+linspace(0,L*b,num=L)
		return(h)
	
	def E(d,b):
		return(b)
	
	def Q(d,b):
		if nanmin(d)==nanmax(d):
			h=linspace(nanmin(d)-0.5,nanmax(d)+0.5,num=b)
		else:
			from math import ceil,floor
			from numpy import array,concatenate,cumsum,ones,sort
			d=sort(d)
			L=len(d)
			w_l=floor(L/b)
			w_h=ceil(L/b)
			n_l=b*ceil(L/b)-L
			n_h=b-n_l
			n=concatenate([ones(ceil(n_l/2))*w_l,ones(n_h)*w_h,ones(floor(n_l/2))*w_l]).astype('int32')
			h=array([nanmin(d)]+[(d[i-1]+d[i])/2 for i in cumsum(n)[:-1]]+[nanmax(d)])
		return(h)
	
	if log:
		from numpy import log10
		data=log10(data)
	if btype is None:
		from numpy import ndarray
		bdict={int:'number',float:'width',ndarray:'edges'}
		btype=bdict[type(bins)]
	bfunc={'number':N,'width':W,'edges':E,'equal':Q}
	hist_bins=bfunc[btype](data,bins)
	plot_bins=hist_bins*1.0
	if log:
		plot_bins=10**plot_bins
	return(data,hist_bins,plot_bins)

def dict_splicer(plot_dict,Ld,Lx):
	"""Dictionary constructor for plotting
	
	This function is a base-level function used by most other plotting functions to construct a list of dictionaries,
	each containing the passed arguments to the underlying plotting calls for each different dataset.
	
	Parameters
	----------
	plot_dict : dict
		Contains the parameters to be passed to the underlying plotting function.
	Ld : int
		Number of plots to be made.
	Lx : list
		Contains the lenght of the data array of each plot to be made.
	
	Returns
	-------
	dict_list : list
		List of dictionaries, one for each plot to be made.
	"""
	dict_list=[]
	dict_keys=plot_dict.keys()
	if 'rasterized' not in dict_keys:
		plot_dict['rasterized']=True
	for i in range(Ld):
		temp_dict={}
		for k in dict_keys:
			try:
				temp=(i for i in plot_dict[k])
			except TypeError:
				temp_dict[k]=plot_dict[k]
			else:
				if type(plot_dict[k]) is str or len(plot_dict[k])==Lx[i]:
					temp_dict[k]=plot_dict[k]
				else:
					temp_dict[k]=plot_dict[k][i]
		dict_list.append(temp_dict)
	return(dict_list)

def percent_finder(data,p):
	"""Level finder for percentage contours
	
	This function is a base-level function used by plots_2d.sigma_cont to define the level which contains the requested
	percentage of the data points.
	
	Parameters
	----------
	data : ndarray
		Describes the n-dimensional position of each data point.
	p : float
		Fraction of the data points to be encircled by the contour.
	
	Returns
	-------
	min_value : float
		Level for the contour.
	"""
	from numpy import cumsum,ravel,sort
	from numpy import min as np_min
	from numpy import sum as np_sum
	
	data_sorted=sort(ravel(data))[::-1]
	data_fraction=cumsum(data_sorted)
	data_fraction/=np_sum(data_sorted)
	try:
		min_value=np_min(data_sorted[data_fraction<p])
	except ValueError:
		min_value=np_min(data_sorted)
	return(min_value)

def plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid_control):
	"""New axis handler
	
	This function is a base-level function used by most other plotting functions to set the current Axes instance
	to ax and returns the old Axes instance to be later reverted.
	
	Parameters
	----------
	xlog : None or bool
		If True, the x-axis scale is set to 'log'.
	ylog : None or bool
		If True, the y-axis scale is set to 'log'.
	xlim : None or array-like
		If given, defines the low and high limits for the x-axis. The first two elements must be int and/or float.
	ylim : None or array-like
		If given, defines the low and high limits for the y-axis. The first two elements must be int and/or float.
	title : None or str
		If given, defines the title of the figure.
	xlabel : None or str
		If given, defines the label of the x-axis.
	ylabel : None or str
		If given, defines the label of the y-axis.
	xinvert : None or bool
		If True, ensures the x-axis is inverted. If False, ensures the x-axis is not inverted.
	yinvert : None or bool
		If True, ensures the y-axis is inverted. If False, ensures the y-axis is not inverted.
	grid_control : None or bool
		If True, ensures the grid is turned on. If False, ensures the grid is turned off.
	
	Returns
	-------
	None
	"""
	from .defaults import Params
	from matplotlib.pyplot import gca,grid,xscale,yscale
	from matplotlib.pyplot import title as plt_title
	from matplotlib.pyplot import xlabel as plt_xlabel
	from matplotlib.pyplot import xlim as plot_xlim
	from matplotlib.pyplot import ylabel as plt_ylabel
	from matplotlib.pyplot import ylim as plt_ylim
	
	if xlog:
		xscale('log')
	if ylog:
		yscale('log')
	if xlim is not None:
		plt_xlim(xlim)
	if ylim is not None:
		plt_ylim(ylim)
	if title is not None:
		plt_title(title)
	if xlabel is not None:
		plt_xlabel(xlabel)
	if ylabel is not None:
		plt_ylabel(ylabel)
	if xinvert:
		if not gca().xaxis_inverted():
			gca().invert_xaxis()
	if yinvert:
		if not gca().yaxis_inverted():
			gca().invert_yaxis()
	if grid_control is None:
		grid_control=Params.grid
	grid(b=grid_control,which=Params.grid_which,axis=Params.grid_axis)

