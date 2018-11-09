def base_hist2D(x,y,c,bin_num,norm,density,cstat,xlog,ylog):
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
		Z=np.histogram2d(x,y,bins=[X,Y],density=density)[0]
		Z*=norm
	return(X,Y,Z)

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
	import matplotlib.pyplot as plot
	
	if xlog:
		plot.xscale('log')
	if ylog:
		plot.yscale('log')
	if xlim is not None:
		plot.xlim(xlim)
	if ylim is not None:
		plot.ylim(ylim)
	if title is not None:
		plot.title(title)
	if xlabel is not None:
		plot.xlabel(xlabel)
	if ylabel is not None:
		plot.ylabel(ylabel)
	if xinvert:
		plot.gca().invert_xaxis()
	if yinvert:
		plot.gca().invert_yaxis()
	plot.grid()

