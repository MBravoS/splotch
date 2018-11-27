#### Definition of all wrappers for 1D plotting

#Histogram
def hist(data,bin_num=None,dens=True,norm=1,c=None,xinvert=False,xlim=None,ylim=None,yinvert=False,xlog=False,ylog=True,
			title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,multi=False):
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from .base_func import plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(data) is not list:
		data=[data]
	L=len(data)
	if bin_num is None:
		bin_num=[int((len(d))**0.4) for d in data]
	if type(bin_num) is not list:
		bin_num=[bin_num+1]*L
	if dens:
		if type(norm)!=int and type(norm)!=float:
			norm=[norm[i]*len(data[i]) for i in range(L)]
		else:
			norm=[norm*len(d) for d in data]
	else:
		if type(norm)!=int and type(norm)!=float:
			norm=[norm[i] for i in range(L)]
		else:
			norm=[norm for d in data]
	if c is None or type(c)==str:
		c=[c]*L
	if type(plabel) is not list:
		plabel=[plabel]*L
	for i in range(L):
		if xlog:
			bins=np.logspace(np.log10(np.nanmin(data[i])),np.log10(np.nanmax(data[i])),num=bin_num[i])
		else:
			if np.nanmin(data[i])==np.nanmax(data[i]):
				bins=np.linspace(np.nanmin(data[i])-0.5,np.nanmax(data[i])+0.5,num=bin_num[i])
			else:
				bins=np.linspace(np.nanmin(data[i]),np.nanmax(data[i]),num=bin_num[i])
		y,x=np.histogram(data[i],bins=bins,density=dens)
		y*=norm[i]
		plt.plot((x[0:-1]+x[1:])/2,y,label=plabel[i],rasterized=True)
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(old_axes)

#Step histogram
def histstep(data,bin_num=None,dens=True,hist_type='step',c='k',xinvert=False,xlim=None,ylim=None,yinvert=False,xlog=False,ylog=True,
			title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,ax=None,multi=False):
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plt
	from .base_func import plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(data) is not list:
		data=[data]
	L=len(data)
	if bin_num is None:
		bin_num=[int((len(d))**0.4) for d in data]
	if type(bin_num) is not list:
		bin_num=[bin_num+1]*L
	if type(c)==str:
		c=[c]*L
	if plabel is None:
		plabel=[plabel]*L
	for i in range(L):
		if xlog:
			bins=np.logspace(np.log10(np.nanmin(data[i])),np.log10(np.nanmax(data[i])),num=bin_num[i])
		else:
			if np.nanmin(data[i])==np.nanmax(data[i]):
				bins=np.linspace(np.nanmin(data[i])-0.5,np.nanmax(data[i])+0.5,num=bin_num[i])
			else:
				bins=np.linspace(np.nanmin(data[i]),np.nanmax(data[i]),num=bin_num[i])
		plt.hist(data[i],bins=bins,density=dens,histtype=hist_type,color=c,label=plabel[i],rasterized=True)
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(old_axes)

# Generalized lines
def line(x,y,n=10,a=1,line_style='solid',c='k',xinvert=False,yinvert=False,cinvert=False,xlog=False,ylog=False,xlim=None,ylim=None,xlabel=None,ylabel=None,plabel=None,
			title=None,lab_loc=0,ax=None,multi=False):
	import numpy as np
	import matplotlib.pyplot as plt
	from .base_func import plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if xlog:
		x=np.logspace(np.log10(x[0]),np.log10(x[1]),num=n)
	else:
		x=np.linspace(x[0],x[1],num=n)
	if ylog:
		y=np.logspace(np.log10(y[0]),np.log10(y[1]),num=n)
	else:
		y=np.linspace(y[0],y[1],num=n)
	plt.plot(x,y,color=c,alpha=a,linestyle=line_style,rasterized=True,label=plabel)
	if plabel is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(ax)

#Plots
def plot(x,y,a=1,line_style='solid',line_colour=None,marker_edge_colour='k',marker_edge_width=0,marker_face_colour='k',marker_size=0,marker_type='o',
			xinvert=False,yinvert=False,cinvert=False,xlog=False,ylog=False,xlim=None,ylim=None,xlabel=None,ylabel=None,plabel=None,
			title=None,lab_loc=0,ax=None,multi=False):
	import numpy as np
	import matplotlib.pyplot as plt
	from .base_func import plot_finalizer
	
	if ax is not None:
		old_axes=axes_handler(ax)
	if type(x) is not list:
		x=[x]
	if type(y) is not list:
		y=[y]
	L=len(x)
	if type(a) is not list:
		a=[a]*L
	if type(line_colour) is not list:
		line_colour=[line_colour]*L
	if type(line_style) is not list:
		line_style=[line_style]*L
	if type(marker_edge_colour) is not list:
		marker_edge_colour=[marker_edge_colour]*L
	if type(marker_edge_width) is not list:
		marker_edge_width=[marker_edge_width]*L
	if type(marker_face_colour) is not list:
		marker_face_colour=[marker_face_colour]*L
	if type(marker_size) is not list:
		marker_size=[marker_size]*L
	if type(marker_type) is not list:
		marker_type=[marker_type]*L
	if plabel is None:
		plabel=[plabel]*L
	for i in range(L):
		plt.plot(x[i],y[i],alpha=a[i],label=plabel[i],linestyle=line_style[i],color=line_colour[i],markeredgecolor=marker_edge_colour[i],
					markeredgewidth=marker_edge_width[i],markerfacecolor=marker_face_colour[i],markersize=marker_size[i],marker=marker_type[i],rasterized=True)
	if plabel[0] is not None:
		plt.legend(loc=lab_loc)
	if not multi:
		plot_finalizer(xlog,ylog,xlim,ylim,title,xlabel,ylabel,xinvert,yinvert)
	if ax is not None:
		old_axes=axes_handler(old_axes)

