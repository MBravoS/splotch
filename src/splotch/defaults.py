########################################################################
########################### SPLOTCH defaults ###########################
########################################################################

####################################
# Default values class
####################################
class Params:
	#Contours
	cont_filled=False
	cont_cmap='viridis'
	cont_linestyle=['solid','dashed','dotted']
	contp_percent=[68.27,95.45]
	contp_output=False
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
