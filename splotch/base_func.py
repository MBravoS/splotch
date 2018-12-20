def axes_handler(new_axis):
	import matplotlib.pyplot as plt
	
	curr_axis=plt.gca()
	plt.sca(new_axis)
	return(curr_axis)

def base_hist2D(x,y,c,bin_num,norm,dens,cstat,xlog,ylog):
	import numpy as np
	import scipy.stats as stats
	
	if xlog:
		X=np.logspace(np.log10(np.nanmin(x)),np.log10(np.nanmax(x)),num=bin_num)
	else:
		if np.nanmin(x)==np.nanmax(x):
			X=np.linspace(np.nanmin(x)-0.5,np.nanmax(x)+0.5,num=bin_num)
		else:
			X=np.linspace(np.nanmin(x),np.nanmax(x),num=bin_num)
	if ylog:
		Y=np.logspace(np.log10(np.nanmin(y)),np.log10(np.nanmax(y)),num=bin_num)
	else:
		if np.nanmin(y)==np.nanmax(y):
			Y=np.linspace(np.nanmin(y)-0.5,np.nanmax(y)+0.5,num=bin_num)
		else:
			Y=np.linspace(np.nanmin(y),np.nanmax(y),num=bin_num)
	if cstat:
		Z=stats.binned_statistic_2d(x,y,c,statistic=cstat,bins=[X,Y])[0]
	else:
		Z=np.histogram2d(x,y,bins=[X,Y],density=dens)[0]
		if dens and norm:
			Z*=1.0*len(x)/norm
	return(X,Y,Z)

def dict_splicer(plot_dict,L):
	dict_list=[]
	dict_keys=plot_dict.keys()
	if 'rasterized' not in dict_keys:
		plot_dict['rasterized']=True
	for i in range(L):
		temp_dict={}
		for k in dict_keys:
			try:
				temp=(i for i in plot_dict[k])
			except TypeError:
				temp_dict[k]=plot_dict[k]
			else:
				if type(plot_dict[k]) is str:
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

