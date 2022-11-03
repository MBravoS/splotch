########################################################################
######### Base functions for axis_funcs, plots_1d and plots_2d #########
########################################################################
from numbers import Number
from math import ceil, floor
from distutils.util import strtobool
from numpy import array, asarray, concatenate, cumsum, empty, histogram2d, linspace, log10, nanmax
from numpy import nanmin, nansum, ndarray, ones, ravel, size, sort
from scipy.stats import binned_statistic_2d
from matplotlib import rcParams
from matplotlib.pyplot import fill_between, fill_betweenx, gca, grid, sca, title as plt_title
from matplotlib.pyplot import xlabel as plt_xlabel, xlim as plt_xlim, xscale
from matplotlib.pyplot import ylabel as plt_ylabel, ylim as plt_ylim, yscale

####################################
# Boolean, unsigned integer, signed integer, float, complex.
####################################
_NUMERIC_KINDS = set('buifc')


####################################
# Generate new axes, return previous
####################################
def axes_handler(new_axis):
    """New axis handler

    Base-level function used by most other plotting functions to set the current Axes instance to
    ax and returns the old Axes instance to be later reverted.

    Parameters
    ----------
    new_axis : Axes object
        The new Axes instance to be set as the current axis.

    Returns
    -------
    curr_axis : Axes object
        The previous Axes instance.
    """

    curr_axis = gca()  # get current axis if it exists, else create one.
    sca(new_axis)
    return(curr_axis)


####################################
# Base function for 2D histograms
####################################
def basehist2D(x, y, c, weights, mbin_type, bin_num, norm, dens, cstat, xlog, ylog):
    """2D histogram base calculation

    Base-level function used by plots_2d.hist2D() and plots_2d.sigma_cont() to calculate the
    underlying histogram.

    Parameters
    ----------
    x : ndarray
        Position of data points in the x axis.
    y : ndarray
        Position of data points in the y axis.
    c : ndarray
        If a valid argument is given in cstat, defines the value used for the binned statistics.
    bin_type : {'number','width','edges','equal'}
        Defines how is understood the value given in bins: 'number' for givinf the desired number
        of bins, 'width' for the width of the bins, 'edges' for the edges of bins, and 'equal' for
        making bins with equal number of elements (or as close as possible). If not given it is
        inferred from the data type of bins: 'number' if int, 'width' if float and 'edges' if ndarray.
    bin_num : int, float or array-like
        Gives the values for the bins, according to bin_type.
    norm : float
        Normalization of the counts.
    dens : bool
        If false the histogram returns raw counts.
    cstat : str or function
        Must be one of the valid str arguments for the statistics variable in scipy.stats.binned_statistic_2d
        ('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or a function that takes a 1D array and
        outputs an integer or float.
    xlog : bool
        If True the scale of the x-axis is logarithmic.
    ylog : bool
        If True the scale of the x-axis is logarithmic.

    Returns
    -------
    x_bins_plot : ndarray
        The bin edges on the x-axis.
    y_bins_plot : ndarray
        The bin edges on the y-axis.
    Z : ndarray
        The value of each bin.
    """

    x_temp, x_bins_hist, x_bins_plot = bin_axis(x, mbin_type[0], bin_num[0], log=xlog)
    y_temp, y_bins_hist, y_bins_plot = bin_axis(y, mbin_type[1], bin_num[1], log=ylog)
    if cstat:
        Z = binned_statistic_2d(x_temp, y_temp, c, statistic=cstat, bins=[x_bins_hist, y_bins_hist])[0]
    else:
        Z = histogram2d(x_temp, y_temp, bins=[x_bins_hist, y_bins_hist], weights=weights,
                        density=False if (size(x_temp) == 0 or size(y_temp) == 0) else dens)[0]
        if norm:
            Z /= norm
            if dens:
                #  bin_area = empty(Z.shape) # commenting out as bin_area not used
                for i in range(len(x_bins_hist) - 1):
                    for j in range(len(y_bins_hist) - 1):
                        Z[i, j] /= (x_bins_hist[i + 1] - x_bins_hist[i]) * (y_bins_hist[j + 1] - y_bins_hist[j])
                scale_temp = (x_temp > min(x_bins_hist)) & (x_temp < max(x_bins_hist))
                scale_temp &= (y_temp > min(y_bins_hist)) & (y_temp < max(y_bins_hist))
                Z *= sum(scale_temp)
    return(x_bins_plot, y_bins_plot, Z)


def bin_axis(data, btype, bins, log=False, plot_centre=False):
    """Bin construction for histograms

    Base-level function used by all histogram-related functions to construct the bins.

    Parameters
    ----------
    data : ndarray
        Data to be binned.
    bin_type : {'number','width','edges','equal'}
        Defines how is understood the value given in bins: 'number' for givinf the desired number
        of bins, 'width' for the width of the bins, 'edges' for the edges of bins, and 'equal' for
        making bins with equal number of elements (or as close as possible). If not given it is
        inferred from the data type of bins: 'number' if int, 'width' if float and 'edges' if ndarray.
    bins : int, float, array-like or list
        Gives the values for the bins, according to bin_type.
    log: bool, optional
        If True, the bins are constructed in logarithmic space and the logarithm of the data is returned.
    step: bool, optional
        If True, returns the edges of the bins, instead of the midpoint values.

    Returns
    -------
    data : ndarray
        The data for the histogram.
    hist_bins : ndarray
        The bin edges.
    plot_bins: ndarray
        The bin centers.

    """

    def N(d, b):
        if (size(d) == 0):  # If no data given.
            return linspace(0.0, 1.0, num=b + 1)
        if nanmin(d) == nanmax(d):
            h = linspace(nanmin(d) - 0.5, nanmax(d) + 0.5, num=b + 1)
        else:
            h = linspace(nanmin(d), nanmax(d), num=b + 1)

        return(h)

    def W(d, b):
        if (size(d) == 0):  # If no data given.
            return linspace(0, 1.0, num=ceil(1.0 / b))
        if nanmin(d) == nanmax(d):
            h = array([nanmin(d) - b / 2, nanmax(d) + b / 2])
        else:
            L = ceil((nanmax(d) - nanmin(d)) / b)
            h = nanmin(d) + linspace(0, L * b, num=L)

        return(h)

    def E(d, b):
        return(b)

    def Q(d, b):
        if not isinstance(b, int):
            raise TypeError('bins must be integer when bin_type="equal"')
        if len(d) < b:
            raise IndexError('the number of bins must be smaller than the length of the data when bin_type="equal"')
        if (size(d) == 0):  # If no data given.
            return linspace(0.0, 1.0, num=b + 1)
        if nanmin(d) == nanmax(d):
            h = linspace(nanmin(d) - 0.5, nanmax(d) + 0.5, num=b + 1)
        else:
            d = sort(d)
            L = len(d)
            w_l = floor(L / b)
            w_h = ceil(L / b)
            n_l = b * ceil(L / b) - L
            n_h = b - n_l
            n = concatenate([ones(ceil(n_l / 2)) * w_l, ones(n_h) * w_h, ones(floor(n_l / 2)) * w_l]).astype('int32')
            h = array([nanmin(d)] + [(d[i - 1] + d[i]) / 2 for i in cumsum(n)[:-1]] + [nanmax(d)])

        return(h)

    if log:
        data = log10(data)

    if btype is None:
        bdict = {int: 'number', float: 'width', ndarray: 'edges'}
        btype = bdict[type(bins)]

    bfunc = {'number': N, 'width': W, 'edges': E, 'equal': Q}
    hist_bins = bfunc[btype](data, bins)
    plot_bins = None if hist_bins is None else hist_bins * 1.0

    if plot_centre:
        plot_bins = (plot_bins[:-1] + plot_bins[1:]) / 2

    if log:
        plot_bins = 10**plot_bins

    return(data, hist_bins, plot_bins)


####################################
# Distribute kwargs dicts
####################################
def dict_splicer(plot_dict, Ld, Lx):
    """Dictionary constructor for plotting

    Base-level function used by most other plotting functions to construct a list of dictionaries,
    each containing the passed arguments to the underlying plotting calls for each different dataset.

    Parameters
    ----------
    plot_dict : dict
        Contains the parameters to be passed to the underlying plotting function.
    Ld : int
        Number of plots to be made.
    Lx : list
        Contains the lenght of the data array of each plot to be made.

    Returns
    -------
    dict_list : list
        List of dictionaries, one for each plot to be made.
    """

    dict_list = []
    dict_keys = plot_dict.keys()
    for i in range(Ld):
        temp_dict = {}
        for k in dict_keys:
            try:
                _ = (i for i in plot_dict[k])
            except TypeError:
                temp_dict[k] = plot_dict[k]
            else:
                if isinstance(plot_dict[k], str) or len(plot_dict[k]) == Lx[i]:
                    temp_dict[k] = plot_dict[k]
                elif (k in 'color' or 'color' in k) and (isinstance(plot_dict[k], str) or isinstance(plot_dict[k][0], Number)):
                    temp_dict[k] = plot_dict[k]
                else:
                    temp_dict[k] = plot_dict[k][i]
        dict_list.append(temp_dict)

    return(dict_list)

####################################
# Density/scaled counts for hexbin
####################################
#def hexarea(value,norm):
#    return(value/norm)
#
#def hexscaled(value,scale):
#    return(1.0*value/scale)


####################################
# General check for numeric values
####################################
def is_numeric(array):
    """Base-level numeric classificator of arrays

    Determine whether the argument has a numeric datatype, when converted to a NumPy array.
    Booleans, unsigned integers, signed integers, floats and complex numbers are the kinds of
    numeric datatype.

    Parameters
    ----------
    array : array-like
        The array to check.

    Returns
    -------
    is_numeric : `bool`
        True if the array has a numeric datatype, False if not.

    Credit for code: StackExchange users Gareth Rees and MSeifert.
    https://codereview.stackexchange.com/questions/128032/check-if-a-numpy-array-contains-numerical-data
    """

    return asarray(array).dtype.kind in _NUMERIC_KINDS


####################################
# General check for numbers
####################################
def is_number(var):
    """Base-level numeric classificator of variables

    Determine whether the argument is a literal number. Booleans, unsigned integers, signed integers,
    floats and complex numbers are considered numbers. Arrays and list-like objects will return False
    even if their data type is numeric in kind.

    Parameters
    ----------
    array : array-like
        The array to check.

    Returns
    -------
    is_numeric : `bool`
        True if the array has a numeric datatype, False if not.

    Credit for code: StackExchange users Gareth Rees and MSeifert.
    https://codereview.stackexchange.com/questions/128032/check-if-a-numpy-array-contains-numerical-data
    """
    from numbers import Number

    return isinstance(var, Number)


####################################
# Levels for contourp
####################################
def percent_finder(data, p):
    """Level finder for percentage contours

    This function is a base-level function used by plots_2d.sigma_cont to define the level which
    contains the requested percentage of the data points.

    Parameters
    ----------
    data : ndarray
        Describes the n-dimensional position of each data point.
    p : float
        Fraction of the data points to be encircled by the contour.

    Returns
    -------
    min_value : float
        Level for the contour.
    """

    data_sorted = sort(ravel(data))[::-1]
    data_fraction = cumsum(data_sorted)
    data_fraction /= nansum(data_sorted)

    try:
        min_value = nanmin(data_sorted[data_fraction < p])
    except ValueError:
        min_value = nanmin(data_sorted)
    return(min_value)


####################################
# Simpler version for curve params
####################################
def simpler_dict_splicer(plot_dict, Ld, Lx):
    """Simpler dictionary constructor specifically for plots_1d.curve() plotting

    Base-level function used to construct a list of dictionaries for the characters in the expressions
    to be replaced with numerical values in plots_1d.curve().

    Parameters
    ----------
    plot_dict : dict
        Contains the parameters to be passed to the underlying plotting function.
    Ld : int
        Number of plots to be made.
    Lx : list
        Contains the length of the data array of each plot to be made.

    Returns
    -------
    dict_list : list
        List of dictionaries, one for each plot to be made.
    """

    dict_list = []
    dict_keys = plot_dict.keys()
    for i in range(Ld):
        temp_dict = {}
        for k in dict_keys:
            try:
                _ = (i for i in plot_dict[k])
            except TypeError:
                temp_dict[k] = plot_dict[k]
            else:
                temp_dict[k] = plot_dict[k][i]
        dict_list.append(temp_dict)
    return(dict_list)


####################################
# Set labels, limits and more
####################################
def plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid_control):
    """New axis handler

    This function is a base-level function used by most other plotting functions to set the current
    Axes instance to ax and returns the old Axes instance to be later reverted.

    Parameters
    ----------
    xlog : None or bool
        If True, the x-axis scale is set to 'log'.
    ylog : None or bool
        If True, the y-axis scale is set to 'log'.
    xlim : None or array-like
        If given, defines the low and high limits for the x-axis. The first two elements must be int
        and/or float.
    ylim : None or array-like
        If given, defines the low and high limits for the y-axis. The first two elements must be int
        and/or float.
    title : None or str
        If given, defines the title of the figure.
    xlabel : None or str
        If given, defines the label of the x-axis.
    ylabel : None or str
        If given, defines the label of the y-axis.
    xinvert : None or bool
        If True, ensures the x-axis is inverted. If False, ensures the x-axis is not inverted.
    yinvert : None or bool
        If True, ensures the y-axis is inverted. If False, ensures the y-axis is not inverted.
    grid_control : None or bool
        If True, ensures the grid is turned on. If False, ensures the grid is turned off.

    Returns
    -------
    None
    """

    if xlog:
        xscale('log')

    if ylog:
        yscale('log')

    if xlim is not None:
        plt_xlim(xlim)
    else:
        plt_xlim(auto=True)

    if ylim is not None:
        plt_ylim(ylim)
    else:
        plt_ylim(auto=True)

    if title is not None:
        plt_title(title)

    if xlabel is not None:
        plt_xlabel(xlabel)

    if ylabel is not None:
        plt_ylabel(ylabel)

    if xinvert:
        if not gca().xaxis_inverted():
            gca().invert_xaxis()

    if yinvert:
        if not gca().yaxis_inverted():
            gca().invert_yaxis()

    if grid_control is not None:
        grid(b=grid_control, which=rcParams['axes.grid.which'], axis=rcParams['axes.grid.axis'])
    else:
        grid(b=rcParams['axes.grid'], which=rcParams['axes.grid.which'], axis=rcParams['axes.grid.axis'])


####################################
# Modified fill_between for hist
####################################
### This function is now not needed as step_filler() has been replaced by stairs() in splt.hist()
def step_filler(x, y, **kwargs):
    """Wrapper for a specific fill_between plot.

    Parameters
    ----------
    x : array-like
        The x-axis values.
    y : array-like
        The y-axis values.
    Other arguments : **kwargs
        Optional kwargs supported by fill_between. Note that it will conflict if 'step' is given.
    Returns
    -------
    val : Recasted value
    """

    temp_y = empty(len(y) + 1)
    temp_y[1:] = y
    temp_y[0] = y[0]
    fill_between(x, temp_y, step='pre', **kwargs)

    return(None)


####################################
# Modified fill_between for hist
####################################
### This function is now not needed as step_fillerx() has been replaced by stairs(..., orientation='horizontal') in splt.hist()
def step_fillerx(x, y, **kwargs):
    """Wrapper for a specific fill_between plot.

    Parameters
    ----------
    x : array-like
        The x-axis values.
    y : array-like
        The y-axis values.
    Other arguments : **kwargs
        Optional kwargs supported by fill_betweenx. Note that it will conflict if 'step' is given.
    Returns
    -------
    val : Recasted value
    """

    temp_y = empty(len(y) + 1)
    temp_y[1:] = y
    temp_y[0] = y[0]
    fill_betweenx(x, temp_y, step='pre', **kwargs)

    return(None)




####################################
# Variable type check
####################################
def val_checker(val):
    """Type check of variable to pass to style functions.

    Parameters
    ----------
    val : str
        String containing the value to evaluate
    Returns
    -------
    val : Recasted value
    """

    if "'" not in val and '"' not in val:
        try:
            val = float(val) if '.' in val else int(val)
        except ValueError:
            try:
                val = bool(strtobool(val))
            except ValueError:
                val = val

    return(val)
