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

def base_hist2D(x,y,c,bin_type,bin_num,norm,dens,cstat,xlog,ylog):
	import numpy as np
	import scipy.stats as stats
	
	x_temp,x_bins_hist,x_bins_plot=binned_axis(x,bin_type[0],bin_num[0],log=xlog)
	y_temp,y_bins_hist,y_bins_plot=binned_axis(y,bin_type[1],bin_num[1],log=ylog)
	if cstat:
		Z=stats.binned_statistic_2d(x_temp,y_temp,c,statistic=cstat,bins=[x_bins_hist,y_bins_hist])[0]
	else:
		Z=np.histogram2d(x_temp,y_temp,bins=[x_bins_hist,y_bins_hist],density=dens)[0]
		if dens and norm:
			Z*=1.0*len(x)/norm
	return(x_bins_plot,y_bins_plot,Z)

def binned_axis(data,btype,bins,log=False):
	import math
	import numpy as np
	
	def N(d,b):
		if np.nanmin(d)==np.nanmax(d):
			h=np.linspace(np.nanmin(d)-0.5,np.nanmax(d)+0.5,num=b)
		else:
			h=np.linspace(np.nanmin(d),np.nanmax(d),num=b+1)
		return(h)
	def W(d,b):
		if np.nanmin(d)==np.nanmax(d):
			h=np.array([np.nanmin(d)-b/2,np.nanmax(d)+b/2])
		else:
			L=math.ceil((np.nanmax(d)-np.nanmin(d))/b)
			h=np.nanmin(d)+np.linspace(0,L*b,num=L)
		return(h)
	def E(d,b):
		return(b)
	def Q(d,b):
		if np.nanmin(d)==np.nanmax(d):
			h=np.linspace(np.nanmin(d)-0.5,np.nanmax(d)+0.5,num=b)
		else:
			d=np.sort(d)
			L=len(d)
			w_l=math.floor(L/b)
			w_h=math.ceil(L/b)
			n_l=b*math.ceil(L/b)-L
			n_h=b-n_l
			n=np.concatenate([np.ones(math.ceil(n_l/2))*w_l,np.ones(n_h)*w_h,np.ones(math.floor(n_l/2))*w_l]).astype('int32')
			h=np.array([np.nanmin(d)]+[(d[i-1]+d[i])/2 for i in np.cumsum(n)[:-1]]+[np.nanmax(d)])
		return(h)
	
	if log:
		data=np.log10(data)
	if btype is None:
		bdict={int:'number',float:'width',np.ndarray:'edges'}
		btype=bdict[type(bins)]
	bfunc={'number':N,'width':W,'edges':E,'equal':Q}
	hist_bins=bfunc[btype](data,bins)
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

def plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert,grid_control):
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
		if not plt.gca().xaxis_inverted():
			plt.gca().invert_xaxis()
	if yinvert:
		if not plt.gca().yaxis_inverted():
			plt.gca().invert_yaxis()
	if grid_control is None:
		from .defaults import Params
		grid_control=Params.grid
	plt.grid(b=grid_control,which='major',axis='both')

