from .plots_1d import *
from .plots_2d import *
import matplotlib.pyplot as plt

### The following code is commented out, it should add the style file but it doesn't ###

## Default styling
#def copy_style():
#  import os
#  import matplotlib
#  from pkg_resources import resource_string
#
#  files=['styles/example.mplstyle',]
#
#  for fname in files:
#    path = os.path.join(matplotlib.get_configdir(),fname)
#    text = resource_string(__name__,fname).decode()
#    open(path,'w').write(text)
#	print(path)
#plot.style.use('astro')
#
## Function to return to default style or upload your own
#def change_style(file=None):
#	import matplotlib as mpl
#	
#	mpl.style.use('default')
#	if file is not None:
#		plot.style.use(file)

