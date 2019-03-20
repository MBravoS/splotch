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
	import matplotlib.pyplot as plt
	
	curr_axis=plt.gca()
	plt.sca(new_axis)
	return(curr_axis)

def base_hist2D(x,y,c,bin_num,norm,dens,cstat,xlog,ylog):
	import numpy as np
	import scipy.stats as stats
	
	x_temp,x_bins_hist,x_bins_plot=binned_axis(x,bin_num=bin_num[0],log=xlog)
	y_temp,y_bins_hist,y_bins_plot=binned_axis(y,bin_num=bin_num[1],log=ylog)
	if cstat:
		Z=stats.binned_statistic_2d(x_temp,y_temp,c,statistic=cstat,bins=[x_bins_hist,y_bins_hist])[0]
	else:
		Z=np.histogram2d(x_temp,y_temp,bins=[x_bins_hist,y_bins_hist],density=dens)[0]
		if dens and norm:
			Z*=1.0*len(x)/norm
	return(x_bins_plot,y_bins_plot,Z)

def binned_axis(data,bin_num,log=False):
	import numpy as np
	
	if log:
		data=np.log10(data)
	if np.nanmin(data)==np.nanmax(data):
		hist_bins=np.linspace(np.nanmin(data)-0.5,np.nanmax(data)+0.5,num=bin_num)
	else:
		hist_bins=np.linspace(np.nanmin(data),np.nanmax(data),num=bin_num)
	plot_bins=hist_bins*1.0
	if log:
		plot_bins=10**plot_bins
	return(data,hist_bins,plot_bins)

def dict_splicer(plot_dict,Ld,Lx):
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
	import numpy as np
	
	data_sorted=np.sort(np.ravel(data))[::-1]
	data_fraction=np.cumsum(data_sorted)
	data_fraction/=np.sum(data_sorted)
	try:
		min_value=np.min(data_sorted[data_fraction<p])
	except ValueError:
		min_value=np.min(data_sorted)
	return(min_value)

def plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert):
	import matplotlib.pyplot as plt
	
	if xlog:
		plt.xscale('log')
	if ylog:
		plt.yscale('log')
	if xlim is not None:
		plt.xlim(xlim)
	if ylim is not None:
		plt.ylim(ylim)
	if title is not None:
		plt.title(title)
	if xlabel is not None:
		plt.xlabel(xlabel)
	if ylabel is not None:
		plt.ylabel(ylabel)
	if xinvert:
		plt.gca().invert_xaxis()
	if yinvert:
		plt.gca().invert_yaxis()
	plt.grid()

