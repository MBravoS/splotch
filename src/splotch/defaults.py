#Object containing the default values for SPLOTCH functions
class Params:
	#Contours
	cont_filled=False
	sigcont_percent=[68.27,95.45]
	sigcont_cmap='viridis'
	sigcont_clim=[0.33,0.67]
	sigcont_linestyle=['solid','dashed','dotted']
	sigcont_output=False
	#Grid
	grid=False
	grid_which='major'
	grid_axis='both'
	#Histograms
	hist1D_output=False
	hist1D_histtype='step'
	hist1D_yaxis_log=False
	hist2D_caxis_log=False
	hist2D_output=False
	#Images
	img_caxis_log=False

class Params0():
	def __init__(self):
		self._params = {}
	def __setitem__(self, key, val):
		self._params[key] = val
	def __getitem__(self, key):
		return self._params[key]

#params = Params0()
#params['test'] = 'yay!'
