#### Definition of all wrappers for 1D plotting

#Histogram
def hist(data,bin_num=None,dens=True,norm=1,c=None,xinvert=False,xlim=None,ylim=None,yinvert=False,xlog=False,ylog=True,
			title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,multi=False):
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plot
	
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
		y,x=np.histogram(data[i],bins=bins,density=dens)
		y*=norm[i]
		plot.plot((x[0:-1]+x[1:])/2,y,label=plabel[i],rasterized=True)
	if plabel[0] is not None:
		plot.legend(loc=lab_loc)
	if not multi:
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

# Generalized lines
def line(x,y,n=10,a=1,s='solid',c='k',xinvert=False,yinvert=False,cinvert=False,xlog=False,ylog=False,xlim=None,ylim=None,xlabel=None,ylabel=None,plabel=None,
			title=None,lab_loc=0,multi=False):
	import numpy as np
	import matplotlib.pyplot as plot
	
	if xlog:
		x=np.logspace(np.log10(x[0]),np.log10(x[1]),num=n)
	else:
		x=np.linspace(x[0],x[1],num=n)
	if ylog:
		y=np.logspace(np.log10(y[0]),np.log10(y[1]),num=n)
	else:
		y=np.linspace(y[0],y[1],num=n)
	plot.plot(x,y,color=c,alpha=a,linestyle=s,rasterized=True,label=plabel)
	if plabel is not None:
		plot.legend(loc=lab_loc)
	if not multi:
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


