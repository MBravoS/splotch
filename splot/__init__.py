def hist(data,bin_num=None,density=True,norm=1,c=None,xinvert=False,xlim=None,ylim=None,yinvert=False,xlog=False,ylog=True,
			title=None,xlabel=None,ylabel=None,plabel=None,lab_loc=0,multi=False):
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plot
	
	if type(data) is not list:
		data=[data]
	L=len(data)
	if bin_num is None:
		bin_num=[int((len(d))**0.5) for d in data]
	if type(bin_num) is not list:
		bin_num=[bin_num+1]*L
	if type(norm)!=int and type(norm)!=float:
		norm=[norm[i]*len(data[i]) for i in range(L)]
	else:
		norm=[norm*len(d) for d in data]
	if c is None or type(c)==str:
		c=[c]*L
	if plabel is None:
		plabel=[plabel]*L
	for i in range(L):
		if xlog:
			bins=np.logspace(np.log10(np.min(data[i])),np.log10(np.max(data[i])),num=bin_num[i])
		else:
			if np.min(data[i])==np.max(data[i]):
				bins=np.linspace(np.min(data[i])-0.5,np.max(data[i])+0.5,num=bin_num[i])
			else:
				bins=np.linspace(np.min(data[i]),np.max(data[i]),num=bin_num[i])
		y,x=np.histogram(data[i],bins=bins,density=density)
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

def hist2D(x,y,bin_num=None,density=True,norm=1,c=None,cstat='mean',xlim=None,ylim=None,clim=[None,None],xinvert=False,yinvert=False,cinvert=False,xlog=False,ylog=False,zlog=True,
			title=None,xlabel=None,ylabel=None,clabel=None,lab_loc=0,multi=False):
	import numpy as np
	import scipy.stats as stats
	import matplotlib.colors as clr
	import matplotlib.pyplot as plot
	
	if bin_num is None:
		bin_num=int((len(x))**0.4)
	else:
		bin_num+=1
	norm=norm*len(x)
	if cinvert:
		cmap='viridis_r'
	else:
		cmap='viridis'
	if xlog:
		X=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),num=bin_num)
	else:
		if np.min(x)==np.max(x):
			X=np.linspace(np.min(x)-0.5,np.max(x)+0.5,num=bin_num)
		else:
			X=np.linspace(np.min(x),np.max(x),num=bin_num)
	if ylog:
		Y=np.logspace(np.log10(np.min(y)),np.log10(np.max(y)),num=bin_num)
	else:
		if np.min(y)==np.max(y):
			Y=np.linspace(np.min(y)-0.5,np.max(y)+0.5,num=bin_num)
		else:
			Y=np.linspace(np.min(y),np.max(y),num=bin_num)
	if c is None:
		Z=np.histogram2d(x,y,bins=[X,Y],normed=density)[0]
		Z*=norm
	else:
		Z=stats.binned_statistic_2d(x,y,c,statistic=cstat,bins=[X,Y])[0]
	if zlog:
		plot.pcolormesh(X,Y,Z.T,norm=clr.LogNorm(),vmin=clim[0],vmax=clim[1],cmap=cmap,rasterized=True)
	else:
		plot.pcolormesh(X,Y,Z.T,vmin=clim[0],vmax=clim[1],cmap=cmap,rasterized=True)
	if clabel is not None:
		cbar=plot.colorbar()
		cbar.set_label(clabel)
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
		if cinvert:
			cbar.ax.invert_yaxis()
		plot.grid()

def img(im,x=None,y=None,bin_num=None,xlim=None,ylim=None,clim=[None,None],xinvert=False,yinvert=False,cinvert=False,clog=False,title=None,
		xlabel=None,ylabel=None,clabel=None,lab_loc=0,multi=False):
	import numpy as np
	import matplotlib.colors as clr
	import matplotlib.pyplot as plot
	
	if cinvert:
		cmap='viridis_r'
	else:
		cmap='viridis'
	if x is None:
		x=np.linspace(np.min(im,axis=0),np.max(im,axis=0),len(im[:,0])+1)
	if y is None:
		y=np.linspace(np.min(im,axis=1),np.max(im,axis=1),len(im[0,:])+1)
	plot.pcolormesh(x,y,im,vmin=clim[0],vmax=clim[1],cmap=cmap,rasterized=True)
	if clabel is not None:
		cbar=plot.colorbar()
		cbar.set_label(clabel)
	if not multi:
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
		if cinvert:
			cbar.ax.invert_yaxis()
		plot.grid()		

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

def scat(x,y,marker_size=20,marker_type='o',a=1,xinvert=False,yinvert=False,cinvert=False,xlog=False,ylog=False,xlim=None,ylim=None,xlabel=None,ylabel=None,c=None,
			clabel=None,plabel=None,title=None,lab_loc=0,multi=False):
	import numpy as np
	import matplotlib.pyplot as plot
	
	if type(x) is not list:
		x=[x]
		y=[y]
	colour=[]
	if c is None:
		for i in range(len(x)):
			colour.append(np.random.uniform(low=0,high=1,size=len(x[i])))
	else:
		if c is not list:
			c=[c]
		if type(c[0]) is not str:
			for i in c:
				colour.append(i)
		else:
			colour=c
	if type(plabel) is not list:
		plabel=[plabel]*len(x)
	if type(marker_size) is not list:
		marker_size=[marker_size]*len(x)
	if type(marker_type) is not list:
		marker_type=[marker_type]*len(x)
	if cinvert:
		cmap='viridis_r'
	else:
		cmap='viridis'
	for i in range(len(x)):
		plot.scatter(x[i],y[i],s=marker_size[i],c=colour[i],label=plabel[i],marker=marker_type[i],edgecolors='none',alpha=a,cmap=cmap,rasterized=True)
	if not multi:
		if xlabel is not None:
			plot.xlabel(xlabel)
		if ylabel is not None:
			plot.ylabel(ylabel)
		if clabel is not None:
			cbar=plot.colorbar(fraction=0.046,pad=0.04)
			cbar.set_label(clabel)
		if plabel[0] is not None:
			plot.legend(loc=lab_loc)
		if xlim is not None:
			plot.xlim(xlim)
		if ylim is not None:
			plot.ylim(ylim)
		if title is not None:
			plot.title(title)
		if xlog:
			plot.xscale('log')
		if ylog:
			plot.yscale('log')
		if xinvert:
			plot.gca().invert_xaxis()
		if yinvert:
			plot.gca().invert_yaxis()
		if cinvert:
			cbar.ax.invert_yaxis()
		plot.grid()

