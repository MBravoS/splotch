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
    #Histograms
    hist1D_output=False
    hist1D_histtype='step'
    hist1D_yaxis_log=False
    hist2D_output=False
    hist2D_caxis_log=False
    #Images
    img_caxis_log=False
    img_cmap='Greys_r'
    grid_color = 'white'
    grid_alpha = 1.0
    grid_ls = '--'
    grid_lw = 1.0
