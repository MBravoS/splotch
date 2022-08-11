########################################################################
############## Definition of all wrappers for 2D plotting ##############
########################################################################

####################################
# Level contours
####################################
def contour(z, x=None, y=None, filled=None, xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False, title=None, xlabel=None,
            ylabel=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):
    """Level contour plotting function.

    This is a wrapper for pyplot.contour() and pyplot.contourf().

    Parameters
    ----------
    z : array-like
        The height values to draw the contours.
    x : array-like, optional
        Position of data points in the x axis.
    y : array-like, optional
        Position of data points in the y axis.
    filled: boolean, optional
        If True, draws filled contours. If not given defaults to the value defined in splotch.Params.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool, optional
        If True, inverts the x-axis.
    yinvert : bool, optional
        If True, inverts the y-axis.
    xlog : bool, optional
        If True, the scale of the x-axis is logarithmic. If not given defaults to the value defined in splotch.Params.
    ylog : bool, optional
        If True, the scale of the x-axis is logarithmic. If not given defaults to the value defined in splotch.Params.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    output : boolean, optional
        If True, returns the edges and values of the underlying histogram plus the levels of the contours.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are QuadContourSet properties.
    **kwargs: QuadContourSet properties, optional
        kwargs are used to specify matplotlib specific properties such as cmap, linewidths, hatches, etc.
        The list of available properties can be found here:
        https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.contour.html

    Returns
    -------
    bin_edges_x : array
        The bin edges for the x axis.
    bin_edges_y : array
        The bin edges for the y axis.
    n : array
        The values of the underlying histogram.
    l : array
        The levels for the contours.
    """

    from matplotlib.pyplot import gca, contour, contourf
    from .base_func import axes_handler, plot_finalizer

    # Set the current axis
    if ax is not None:
        old_axes = axes_handler(ax)
    else:
        ax = gca()
        old_axes = ax

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par={**plot_kw, **kwargs} # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    if filled:
        contourset = contourf(x, y, z, **plot_par)
    else:
        contourset = contour(x, y, z, **plot_par)

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)

    return contourset


####################################
# Contours from density histograms
####################################
def contourp(x, y, percent=None, filled=None, bin_type=None, bins=None, smooth=0.0, max_spacing=True, xlim=None, ylim=None, xinvert=False, yinvert=False,
             xlog=False, ylog=False, title=None, labels=None, xlabel=None, ylabel=None, lab_loc=0, ax=None, grid=None, output=None, plot_kw={}, **kwargs):
    """Contour function, encircling the highest density regions that contain the given percentages of the sample.

    Parameters
    ----------
    x : array-like
        Position of data points in the x axis.
    y : array-like
        Position of data points in the y axis.
    percent : float or array-like, optional.
        The percentages of the sample for the contours to encircle.
    bin_type : {'number', 'width', 'edges', 'equal'}, optional
        Defines the format for the value(s) given in `bins`.
            * 'number' specifies the desired number of bins.
            * 'width' specifies the width of the bins.
            * 'edges' for the edges of bins
            * 'equal' results in bins with an equal number of elements (or as close as possible).
        If `bin_type` not given, it is inferred from the data type of bins, i.e., 'number' if int,
        'width' if float and 'edges' if a list-like object.
    bins : int, float, array-like, optional
        Specifies the value(s) for the bins, according to bin_type.
    smooth : float, optional
        The standard deviation for the Gaussian kernel. Default: 0.0 (No smoothing).
    max_spacing : boolean, optional
        If True, maximises the separation between colours drawn from the colour map. Default: True.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool, optional
        If True, inverts the x-axis.
    yinvert : bool, optional
        If True, inverts the y-axis.
    xlog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    labels : array-like or boolean, optional
        Specifies label(s) for the contour(s). If given as an array of strings, sets the labels for each contour.
        Must be of equal length to number specified by 'percent'. If False, no legend is drawn.
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    output : boolean, optional
        If True, returns the edges and values of the underlying histogram in addition to the default output of.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are QuadContourSet properties.
    **kwargs: QuadContourSet properties, optional
        kwargs are used to specify matplotlib specific properties such as cmap, linewidths, hatches, etc.
        The list of available properties can be found here:
        https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.contour.html

    Returns
    -------
    bin_edges_x : array
        The bin edges for the x axis.
    bin_edges_y : array
        The bin edges for the y axis.
    n : array
        The values of the underlying histogram.
    l : array
        The levels for the contours.
    """

    from warnings import warn
    from matplotlib import lines, patches, rcParams
    from matplotlib.cm import get_cmap
    from matplotlib.pyplot import gca, contour, contourf, legend  # , Normalize
    from numpy import array, linspace, round, ndarray
    from scipy.ndimage.filters import gaussian_filter

    from .base_func import axes_handler, basehist2D, percent_finder, plot_finalizer
    from .defaults import Params

    # Check deprecated parameters
    if 'plabel' in kwargs.keys():
        warn("'plabel' will be deprecated, use instead 'labels'", DeprecationWarning, stacklevel=2)
        labels = kwargs.pop('plabel')

    # Initialise defaults
    if filled is None:
        filled = Params.cont_filled
    if percent is None:
        percent = Params.contp_percent
    if 'cmap' not in kwargs.keys() and 'cmap' not in plot_kw.keys() and 'colors' not in kwargs.keys() and 'colors' not in plot_kw.keys():
        plot_kw['cmap'] = Params.cont_cmap
    if output is None:
        output = Params.contp_output

    # Set the current axis
    if ax is not None:
        old_axes = axes_handler(ax)
    else:
        ax = gca()
        old_axes = ax

    if not isinstance(percent, ndarray):
        percent = array([percent]).flatten()

    if not isinstance(bin_type, (list, tuple, ndarray)):
        bin_type = [bin_type] * 2

    if not isinstance(bins, (list, tuple)):
        if bins is None:
            bins = max([10, int(len(x)**0.4)])  # Defaults to min of 10 bins
        bins = [bins] * 2

    percent = percent[::-1]  # reverse the order of percent

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par = {**plot_kw, **kwargs} # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    if filled:
        plot_par['extend'] = 'max'

    # if drawing <4 lines with no color specified, use 'solid', 'dashed' and then 'dotted'
    if not filled and len(percent) < 4 and 'colors' not in plot_par.keys():
        plot_par['colors'] = [next(ax._get_lines.prop_cycler)['color']] * len(percent)

    if 'colors' in plot_par.keys():
        if isinstance(plot_par['colors'], str):
            plot_par['colors'] = [plot_par['colors'] for i in range(len(percent))]

        if 'cmap' in plot_par.keys():
            plot_par.pop('cmap')
    elif max_spacing:
        if isinstance(plot_par['cmap'], str):
            plot_par['cmap'] = get_cmap(plot_par['cmap'])

        plot_par['colors'] = [plot_par['cmap'](i) for i in linspace(0, 1, len(percent))]
        plot_par.pop('cmap')

    # if drawing <4 lines with no color specified, use 'solid', 'dashed' and then 'dotted'
    if not filled and len(percent) < 4 and 'linestyles' not in plot_par.keys():
        plot_par['linestyles'] = [['solid', 'dashed', 'dotted'][i] for i in range(len(percent))][::-1]

    # Validate labels array
    if isinstance(labels, (list, tuple, ndarray)):
        if (len(labels) != len(percent)):
            raise ValueError(f"Length of labels ({len(labels)}) does not match length of percent ({len(percent)}).")
    else:
        if labels is None:
            if rcParams['text.usetex']:
                labels = [f'{round(p,1)}\\%' for p in percent]
            else:
                labels = [f'{round(p,1)}%' for p in percent]

    X, Y, Z = basehist2D(x, y, c=None, bin_type=bin_type, bin_num=bins, norm=None, dens=None, cstat=None, xlog=xlog, ylog=ylog)
    X = (X[:-1] + X[1:]) / 2
    Y = (Y[:-1] + Y[1:]) / 2

    level = array([percent_finder(Z, p / 100) for p in percent])
    if filled:  # Removed `func_dict` for readability
        contourset = contourf(X, Y, gaussian_filter(Z.T, sigma=smooth), levels=level, **plot_par)
    else:
        contourset = contour(X, Y, gaussian_filter(Z.T, sigma=smooth), levels=level, **plot_par)

    if labels:
        plot_par['colors'] = contourset.colors
        if isinstance(plot_par['colors'], str):
            plot_par['colors'] = [plot_par['colors']] * len(percent)

        if contourset.linestyles is None:
            plot_par['linestyles'] = ['solid'] * len(percent)
        elif isinstance(contourset.linestyles, str):
            plot_par['linestyles'] = [contourset.linestyles] * len(percent)
        else:
            plot_par['linestyles'] = contourset.linestyles

        if contourset.alpha is None:
            plot_par['alpha'] = [1.0 for i in range(len(percent))]
        elif isinstance(contourset.alpha, (float, int)):
            plot_par['alpha'] = [contourset.alpha for i in range(len(percent))]
        else:
            plot_par['alpha'] = contourset.alpha

        if filled:
            legend([patches.Patch(color=plot_par['colors'][i], alpha=plot_par['alpha'][i]) for i in range(len(percent))],
                   labels, numpoints=1, loc=lab_loc)
        else:
            legend([lines.Line2D([0, 1], [0, 1], color=plot_par['colors'][i], linestyle=plot_par['linestyles'][i], alpha=plot_par['alpha'][i]) for i in range(len(percent))],
                   labels, numpoints=1, loc=lab_loc)

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)

    if output:
        return(contourset, X, Y, Z.T)
    else:
        return(contourset)


####################################
##########  Error bands  ###########
####################################
def errorband(x, y, yerr, line=False, xlim=None, ylim=None,
              xinvert=False, yinvert=False, xlog=False, ylog=None, title=None, xlabel=None, ylabel=None,
              label=None, lab_loc=0, ax=None, grid=None, line_kw={}, band_kw={}, **kwargs):
    """Error line and band plotting function.

    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like.
    y : array-like or list
        If list it is assumed that each elemement is array-like.
    yerr : array-like or list, optional
        Defines the length of the errobars in the y-axis. If list it is assumed that each elemement is array-like.
    line : boolean, optional
        If True, draw a line that follows the statistic defined in line_stat.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If True, inverts the x-axis.
    yinvert : bool or list, optional
        If True, inverts the y-axis.
    xlog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the label for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, antialiasing, etc.
        A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    None
    """

    from splotch.base_func import axes_handler,bin_axis,plot_finalizer

    import numpy as np
    from numbers import Number
    import scipy.stats as stats
    from numpy import array,percentile
    from functools import partial
    import matplotlib.colors as clr
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import gca
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if ylog is None:
        from splotch.defaults import Params
        ylog = Params.hist1D_yaxis_log
    if 'linewidth' not in band_kw.keys():
        band_kw['linewidth'] = 0
    if 'alpha' not in band_kw.keys():
        band_kw['alpha'] = 0.4

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # band_par={**plot_kw, **kwargs} # For Python > 3.5
    band_kw.update(kwargs)

    if len(array(yerr).shape) == 2:
        plt.fill_between(x, y - yerr[0], y + yerr[1], label=label, **band_kw)
    else:
        plt.fill_between(x, y - yerr, y + yerr, label=label, **band_kw)

    if line:
        plt.plot(x, y, **line_kw)
    if label is not None:
        plt.legend(loc=lab_loc)

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)

####################################
# Error bars
####################################
def errorbar(x, y, xerr=None, yerr=None, xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False,
             title=None, xlabel=None, ylabel=None, label=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):
    """Errorbar plotting function.

    This is a wrapper for pyplot.errorbar().

    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like.
    y : array-like or list
        If list it is assumed that each elemement is array-like.
    xerr : array-like or list, optional
        Defines the length of the errobars in the x-axis. If list it is assumed that each elemement is array-like.
    yerr : array-like or list, optional
        Defines the length of the errobars in the y-axis. If list it is assumed that each elemement is array-like.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If True, inverts the x-axis.
    yinvert : bool or list, optional
        If True, inverts the y-axis.
    xlog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the label for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, antialiasing, etc.
        A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    None
    """

    from .base_func import axes_handler, dict_splicer, plot_finalizer

    from matplotlib.pyplot import errorbar, legend, gca
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    # Convert data to lists if needed
    if not isinstance(x, list): x = [x]
    if not isinstance(y, list): y = [y]
    if not isinstance(xerr, list): xerr = [xerr]
    if not isinstance(yerr, list): yerr = [yerr]

    L = len(x)
    if not isinstance(label, list):
        label = [label for i in range(L)]

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par={**plot_kw, **kwargs} # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    for i in range(L):
        errorbar(x[i], y[i], xerr=xerr[i], yerr=yerr[i], label=label[i], **plot_par[i])
    if any(label):
        legend(loc=lab_loc)
    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)

####################################
###########  Error boxes  ##########
####################################
def errorbox(x, y, xerr=None, yerr=None, xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False, box_type='ellipse',
             title=None, xlabel=None, ylabel=None, label=None, grid=None, lab_loc=0, ax=None, plot_kw={}, **kwargs):
    """Errorbox plotting function.

    This is a wrapper around matplotlib PatchCollections with a matplotlib errorbar functionality.

    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like.
    y : array-like or list
        If list it is assumed that each elemement is array-like.
    xerr : array-like or list, optional
        Defines the length of the errobars in the x-axis. If list it is assumed that each elemement is array-like.
    yerr : array-like or list, optional
        Defines the length of the errobars in the y-axis. If list it is assumed that each elemement is array-like.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If True, inverts the x-axis.
    yinvert : bool or list, optional
        If True, inverts the y-axis.
    xlog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    box_type : str
        The type of box to plot, patch types include: ellipse | rectangle (Default: ellipse).
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the label for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Patches properties.
    **kwargs: Patch properties, optional
        kwargs are used to specify matplotlib specific properties such as facecolor, linestyle, alpha, etc.
        A list of available `Patch` properties can be found here:
        https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.Rectangle.html

    Returns
    -------
    None
    """

    from .base_func import axes_handler, dict_splicer, plot_finalizer

    from matplotlib.pyplot import errorbar, legend, gca
    from numpy import shape, full, array
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Ellipse, Rectangle, Patch
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if not isinstance(x, list): x = [x]
    if not isinstance(y, list): y = [y]
    if not isinstance(xerr, list): xerr = [xerr]
    if not isinstance(yerr, list): yerr = [yerr]

    boxdict = {'rec': Rectangle, 'ell': Ellipse}
    if (box_type.lower()[:3] not in ['rec', 'ell']):
        raise ValueError(f"box_type '{box_type}' not recognised.")

    L = len(x)
    if not isinstance(label, list):
        label = [label for i in range(L)]

    # Validate format of xerr and yerr
    for i in range(L):
        # x-axis errors
        if (shape(xerr[i]) == ()):  # Single error for all points
            xerr[i] = full((2, len(x[i])), xerr[i])
        else:
            if (len(shape(xerr[i])) == 1):
                if (shape(xerr[i])[0] == len(x[i])):  # single error for each point
                    xerr[i] = array([xerr[i], xerr[i]])
                elif (shape(xerr[i])[0] == 2):  # separate upper and lower errors for all points
                    xerr[i] = full((len(x[i]), 2), xerr[i]).T
                else:
                    raise ValueError(f"Invalid shape ({shape(xerr[i])}) for 'xerr' array.")
            elif (len(shape(xerr[i])) == 2):  # separate upper and lower errors for each point
                xerr[i] = array(xerr[i])
                if (shape(xerr[i])[0] != 2 or shape(xerr[i])[1] != len(x[i])):
                    raise ValueError(f"Invalid shape ({shape(xerr[i])}) for 'xerr' array.")

        # y-axis errors
        if (shape(yerr[i]) == ()):  # single error for all points
            yerr[i] = full((2, len(y[i])), yerr[i])
        else:
            if (len(shape(yerr[i])) == 1):
                if (shape(yerr[i])[0] == len(y[i])):  # single error for each point
                    yerr[i] = array([yerr[i], yerr[i]])
                elif (shape(yerr[i])[0] == 2):  # separate upper and lower errors for all points
                    yerr[i] = full((len(y[i]), 2), yerr[i]).T
                else:
                    raise ValueError(f"Invalid shape ({shape(yerr[i])}) for 'yerr' array.")
            elif (len(shape(yerr[i])) == 2):  # separate upper and lower errors for each point
                yerr[i] = array(yerr[i])
                if (shape(yerr[i])[0] != 2 or shape(yerr[i])[1] != len(y[i])):
                    raise ValueError(f"Invalid shape ({shape(yerr[i])}) for 'yerr' array.")

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par={**plot_kw, **kwargs} # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    # Loop over data points; create box/ellipse from errors at each point
    boxdict = {'rec': Rectangle, 'ell': Ellipse}
    boxhandles = []
    for i in range(L):
        errorboxes = []
        for j, (xx, yy, xe, ye) in enumerate(zip(x[i], y[i], xerr[i].T, yerr[i].T)):
            errorboxes.append(boxdict[box_type.lower()[:3]]((xx - xe[0], yy - ye[0]), xe.sum(), ye.sum()))

        # Create and add patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes, **plot_par[i])
        boxhandles.append(Patch(**plot_par[i]))
        ax.add_collection(pc)

    if any(label):
        legend(handles=boxhandles, labels=label, loc=lab_loc)
    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)

####################################
# Hexagonal 2D histogram
####################################
def hexbin(x, y, bins=None, binlims=None, dens=True, scale=None,
           c=None, cstat=None, xlim=None, ylim=None, clim=[None, None], nmin=0,
           xinvert=False, yinvert=False, cbar_invert=False, xlog=False, ylog=False, clog=None, title=None, xlabel=None,
           ylabel=None, clabel=None, lab_loc=0, ax=None, grid=None, output=None, plot_kw={}, **kwargs):
    """Hexagonal 2D bins function.

    Parameters
    ----------
    x : array-like
        Position of data points in the x axis.
    y : array-like
        Position of data points in the y axis.
    bins : int or list, optional
        Gives the number of bins
    binlims : array-like, optional
        Defines the limits for the bin range, given as a one-dimensional list containing four elements (left, right, bottom, top).
        Note, if xlog=True or ylog=True, then the values of binlim must also be logged. The default assigns the limits based on
        gridsize, x, y, xscale and yscale.
    dens : bool or list, optional
        If false the histogram returns raw counts.
    c : array-like, optional
        If a valid argument is given in cstat, defines the value used for the binned statistics.
    cstat : str or function, optional
        Must be one of the valid str arguments for the statistics variable in scipy.stats.binned_statistic_2d
        ('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or a function that takes a 1D array and
        outputs an integer or float.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    clim : list, optional
        Defines the limits of the colour map ranges, it must contain two elements (lower and higer limits).
    nmin : int, optional (default: 0)
        The minimum number of points required in a bin in order to be plotted.
    xinvert : bool, optional
        If True, inverts the x-axis.
    yinvert : bool, optional
        If True, inverts the y-axis.
    cbar_invert : bool, optional
        If True, inverts the direction of the colour bar (not the colour map).
    xlog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    clog : bool, optional
        If True, the colour map is changed from linear to logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    clabel : str, optional
        Setting `clabel` triggers the generation of a colourbar with axis label given by its value.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    output : boolean, optional
        If True, returns the edges and values of the histogram.
    plot_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib hexbin function.
        Parameters will be overwritten if also given implicitly as a **kwarg.
    **kwargs : pcolormesh properties, optional
        kwargs are used to specify matplotlib specific properties such as cmap, norm, edgecolors etc.
        https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hexbin.html

    Returns
    -------
    n : array
        The values of the histogram. Only provided if output is True.
    x_edges : array
        The bin edges for the x axis. Only provided if output is True.
    y_edges : array
        The bin edges for the y axis. Only provided if output is True.
    """
    
    from numpy import diff, log10, nan, nanmin, nanmedian, nanmax, nanstd, unique, size, zeros, shape
    from matplotlib.colors import LogNorm
    from matplotlib.pyplot import gca, hexbin, colorbar
    from .base_func import axes_handler, plot_finalizer
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if type(bins) not in [list, tuple]:
        if bins is None:
            bins = max([10, int(len(x)**0.4)])  # Defaults to min of 10 bins
        bins = [bins, int(bins / (3**0.5))]

    if None in (clog, output):
        from .defaults import Params
        if clog is None:
            clog = Params.hist2D_caxis_log
        if output is None:
            output = Params.hist2D_output

    if size([x, y]) == 0:  # Zero-sized arrays given
        if (clog is True): raise ValueError("Cannot set 'clog'=True if zero-size array given.")
        if (cstat is not None): raise ValueError(f"Cannot compute statistic (cstat='{cstat}') on zero-size array, set cstat=None if no data given.")

    # Create a temporary duplicate of x and y data
    temp_x = x * 1.0  # this should probably use: copy.deepcopy(x)
    temp_y = y * 1.0

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par = {**plot_kw, **kwargs}  # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)
    if binlims:
        plot_par['extent'] = binlims

    if c is not None:
        plot_par['C'] = c
    if cstat:
        cstat_func = {'min': nanmin, 'mean': nanmax, 'median': nanmedian, 'max': nanmax, 'std': nanstd}
        if cstat in cstat_func.keys():
            plot_par['reduce_C_function'] = cstat_func[cstat]
        else:
            plot_par['reduce_C_function'] = cstat

    if clim:
        plot_par['vmin'] = clim[0]
        plot_par['vmax'] = clim[1]

    if nmin:
        plot_par['mincnt'] = nmin

    if xlog:
        plot_par['xscale'] = 'log'
        temp_x = log10(temp_x)

    if ylog:
        plot_par['yscale'] = 'log'
        temp_y = log10(temp_y)

    if clog:
        plot_par['bins'] = 'log'
        if 'mincnt' not in plot_par.keys():
            plot_par['mincnt'] = 1

    if dens and c is None:
        # This is nasty, but seems to be the quickest way to do this without fully rewriting hexbin here
        hist_return = hexbin(temp_x, temp_y, gridsize=bins)
        hist_return.remove()
        offsets = hist_return.get_offsets()
        offsets_x = unique(offsets[:, 0])
        offsets_y = unique(offsets[:, 1])
        hex_area = diff(offsets_x)[0] * 2 * diff(offsets_y)[0]

        def density_scaling(bin_data):
            bin_dens = 1.0 * len(bin_data) / (hex_area)
            if scale:
                bin_dens /= 1.0 * scale
            else:
                bin_dens /= len(x)
            return(bin_dens)

        plot_par['C'] = y
        plot_par['reduce_C_function'] = density_scaling

    hist_return = hexbin(x, y, gridsize=bins, **plot_par)

    if clabel is not None:
        cbar = colorbar()
        cbar.set_label(clabel)
        if cbar_invert:
            cbar.ax.invert_yaxis()

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)

    if old_axes is not ax:
        old_axes = axes_handler(old_axes)
    if output:
        return(hist_return.get_array(), hist_return.get_offsets())

########################################
## 2D histogram and binned statistics ##
########################################
def hist2D(x, y, bin_type=None, bins=None, dens=True, scale=None, c=None, cstat=None, xlim=None, ylim=None, clim=[None, None], nmin=0,
           xinvert=False, yinvert=False, cbar_invert=False, xlog=False, ylog=False, clog=None, title=None, xlabel=None,
           ylabel=None, clabel=None, lab_loc=0, ax=None, grid=None, output=None, plot_kw={}, **kwargs):
    """2D histogram function.

    Parameters
    ----------
    x : array-like
        Position of data points in the x axis.
    y : array-like
        Position of data points in the y axis.
    bin_type : {'number','width','edges','equal'}, optional
        Defines how is understood the value given in bins: 'number' for givinf the desired number of
        bins, 'width' for the width of the bins, 'edges' for the edges of bins, and 'equal' for
        making bins with equal number of elements (or as close as possible). If not given it is
        inferred from the data type of bins: 'number' if int, 'width' if float and 'edges' if ndarray.
    bins : int, float, array-like or list, optional
        Gives the values for the bins, according to bin_type.
    dens : bool or list, optional
        If false the histogram returns raw counts.
    scale : float or list, optional
        Scaling of the data counts.
    c : array-like, optional
        If a valid argument is given in cstat, defines the value used for the binned statistics.
    cstat : str or function, optional
        Must be one of the valid str arguments for the statistics variable in scipy.stats.binned_statistic_2d
        ('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or a function that takes a 1D array and
        outputs an integer or float.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    clim : list, optional
        Defines the limits of the colour map ranges, it must contain two elements (lower and higer limits).
    nmin : int, optional (default: 0)
        The minimum number of points required in a bin in order to be plotted.
    xinvert : bool, optional
        If True, inverts the x-axis.
    yinvert : bool, optional
        If True, inverts the y-axis.
    cbar_invert : bool, optional
        If True, inverts the direction of the colour bar (not the colour map).
    xlog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    clog : bool, optional
        If True, the colour map is changed from linear to logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    clabel : str, optional
        Setting `clabel` triggers the generation of a colourbar with axis label given by its value.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    output : boolean, optional
        If True, returns the edges and values of the histogram.
    plot_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib pcolormesh function.
        Parameters will be overwritten if also given implicitly as a **kwarg.
    **kwargs : pcolormesh properties, optional
        kwargs are used to specify matplotlib specific properties such as cmap, norm, edgecolors etc.
        https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html

    Returns
    -------
    n : array
        The values of the histogram. Only provided if output is True.
    x_edges : array
        The bin edges for the x axis. Only provided if output is True.
    y_edges : array
        The bin edges for the y axis. Only provided if output is True.
    """

    from numpy import nan, size, zeros, shape
    from matplotlib.colors import LogNorm
    from matplotlib.pyplot import gca, pcolormesh, colorbar
    from .base_func import axes_handler, basehist2D, plot_finalizer
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if not isinstance(bin_type, list):
        bin_type = [bin_type] * 2
    if type(bins) not in [list, tuple]:
        if bins is None:
            bins = max([10, int(len(x)**0.4)])  # Defaults to min of 10 bins
        bins = [bins] * 2

    if None in (clog, output):
        from .defaults import Params
        if clog is None:
            clog = Params.hist2D_caxis_log
        if output is None:
            output = Params.hist2D_output

    if size([x, y]) == 0:  # Zero-sized arrays given
        if (clog is True): raise ValueError("Cannot set 'clog'=True if zero-size array given.")
        if (cstat is not None): raise ValueError(f"Cannot compute statistic (cstat='{cstat}') on zero-size array, set cstat=None if no data given.")

    X, Y, Z = basehist2D(x, y, c, bin_type, bins, scale, dens, cstat, xlog, ylog)

    # Also get counts for number threshold cut
    if (size([x, y]) == 0):
        counts = zeros(shape=shape(Z))
    else:
        _, _, counts = basehist2D(x, y, c, bin_type, bins, None, False, None, xlog, ylog)

    # Cut bins which do not meet the number count threshold
    Z[counts < nmin] = nan

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par = {**plot_kw, **kwargs}  # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    if clog:
        pcolormesh(X, Y, Z.T, norm=LogNorm(vmin=clim[0], vmax=clim[1], clip=False), **plot_par)
    else:
        pcolormesh(X, Y, Z.T, vmin=clim[0], vmax=clim[1], **plot_par)

    if clabel is not None:
        cbar = colorbar()
        cbar.set_label(clabel)
        if cbar_invert:
            cbar.ax.invert_yaxis()
    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)
    if output:
        return(Z.T, X, Y)


####################################
######  Image from 2D array  #######
####################################
def img(im, x=None, y=None, xlim=None, ylim=None, clim=[None, None], cmin=0, xinvert=False, yinvert=False, cbar_invert=False, clog=None, 
        title=None, xlabel=None, ylabel=None, clabel=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):

    """2D pixel-based image plotting function.

    Parameters
    ----------
    im : array-like
        Value for each pixel in an x-y 2D array, where the first dimension is the x-position and the
        second is the y-position.
    x : array-like, optional
        Position of data points in the x axis.
    y : array-like, optional
        Position of data points in the y axis.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    clim : list, optional
        Defines the limits of the colour map ranges, it must contain two elements (lower and higer limits).
    clog : bool, optional
        If True, the colour map is changed from linear to logarithmic.
    xinvert : bool, optional
        If True, inverts the x-axis.
    yinvert : bool, optional
        If True, inverts the y-axis.
    cbar_invert : bool, optional
        If True, inverts the direction of the colour bar (not the colour map).
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    clabel : str, optional
        Setting `clabel` triggers the generation of a colourbar with axis label given by its value.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib pcolormesh function.
        Parameters will be overwritten if also given implicitly as a **kwarg.
    **kwargs : pcolormesh properties, optional
        kwargs are used to specify matplotlib specific properties such as `cmap`, `marker`, `norm`, etc.
        A list of available `pcolormesh` properties can be found here:
        https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html

    Returns
    -------
    None
    """

    from numpy import arange, meshgrid
    from matplotlib.colors import LogNorm
    from matplotlib.pyplot import gca, pcolormesh, colorbar
    from .base_func import axes_handler, plot_finalizer
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if x is None:
        x = arange(len(im[:, 0]) + 1)
    if y is None:
        y = arange(len(im[0, :]) + 1)
    if clog is None:
        from .defaults import Params
        clog = Params.img_caxis_log

    X, Y = meshgrid(x, y)

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par = {**plot_kw, **kwargs}  # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    if clog:
        pcolormesh(X, Y, im.T, norm=LogNorm(vmin=clim[0], vmax=clim[1], clip=True), **plot_par)
    else:
        pcolormesh(X, Y, im.T, vmin=clim[0], vmax=clim[1], **plot_par)

    if clabel is not None:
        cbar = colorbar()
        cbar.set_label(clabel)
        if cbar_invert:
            cbar.ax.invert_yaxis()

    plot_finalizer(False, False, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)

    if old_axes is not ax:
        old_axes = axes_handler(old_axes)


####################################
##########  Scatter plots ##########
####################################
def scatter(x, y, c=None, xlim=None, ylim=None, clim=None, density=False, xinvert=False, yinvert=False, cbar_invert=False, xlog=False, ylog=False, title=None,
            xlabel=None, ylabel=None, clabel=None, label=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):
    """2D pixel-based image plotting function.

    Parameters
    ----------
    x : array-like or list
        Position of data points in the x-axis.
    y : array-like or list
        Position of data points in the y-axis.
    c : array-like or list or str, optional
        Value of data points in the z-axis (colour-axis).
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    clim : tuple-like, optional
        Defines the limits of the colour-axis, it must contain two elements (lower and higer limits).
        Functions equivalently to the `vmin, vmax` arguments used by `colors.Normalize`. If both are
        given, `clim` takes priority.
    density : bool, optional
        If True, color-codes points by their spatial density to nearby points using a Gaussian
        kernel density estimate. If 'c' also given, 'density' takes precedence. Default: False.
    xinvert : bool, optional
        If True, inverts the x-axis.
    yinvert : bool, optional
        If True, inverts the y-axis.
    cbar_invert : bool, optional
        If True, inverts the direction of the colour bar (not the colour map).
    xlog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool, optional
        If True, the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    clabel : str, optional
        Setting `clabel` triggers the generation of a colourbar with axis label given by its value.
    label : str, optional
        Sets the label for the scatter plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib scatter function.
        Parameters will be overwritten if also given implicitly as a **kwarg.
    **kwargs : Collection properties, optional
        kwargs are used to specify matplotlib specific properties such as cmap, marker, norm, etc.
        A list of available `Collection` properties can be found here:
        https://matplotlib.org/api/collections_api.html#matplotlib.collections.Collection

    Returns
    -------
    paths
        A list of PathCollection objects representing the plotted data.
    """

    from numpy import array, dtype, shape, vstack
    from matplotlib.pyplot import gca, scatter, legend
    from scipy.stats import gaussian_kde
    from warnings import warn

    from .base_func import axes_handler, dict_splicer, plot_finalizer
    from .axis_func import colorbar
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if (not isinstance(x, list)) or (len(shape(x)) == 1 and array(x).dtype is not dtype('O')):
        x = [x]
    if (not isinstance(y, list)) or (len(shape(y)) == 1 and array(y).dtype is not dtype('O')):
        y = [y]

    L = len(x)

    if (not isinstance(c, list)) or (len(shape(c)) == 1 and array(c).dtype is not dtype('O')):
        c = [c]
        if isinstance(c[0], str) or c[0] is None:
            c = [c[0] for i in range(L)]
    if not isinstance(label, list):
        label = [label for i in range(L)]

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par = {**plot_kw, **kwargs}  # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    # Insert clim as vmin, vmax into **kwargs dictionary, if given.
    if (clim is not None):
        try:
            _ = (e for e in clim)
            if (len(clim) == 2):
                plot_par['vmin'] = clim[0]
                plot_par['vmax'] = clim[1]
            else:
                raise TypeError("`clim` must be of iterable type and have two values only.")
        except (TypeError):
            raise TypeError("`clim` must be of iterable type and have two values only.")

    if (density is True):
        if (all([kk is not None for kk in c])):
            warn("Cannot specify both `c` and `density`, ignoring `c`.")
        c = [None] * L

        for i in range(L):
            xy = vstack([x[i], y[i]])
            c[i] = gaussian_kde(xy)(xy)  # Calculate the Gaussian kernel density estimate

    # Create 'L' number of plot kwarg dictionaries to parse into each scatter call
    plot_par = dict_splicer(plot_par, L, [len(i) for i in x])

    paths = []
    for i in range(L):
        p = scatter(x[i], y[i], c=c[i], label=label[i], **plot_par[i])
        paths.append(p)
    if clabel is not None:
        cbar = colorbar()
        cbar.set_label(clabel)
        if cbar_invert:
            cbar.ax.invert_yaxis()
    if any(label):
        legend(loc=lab_loc)

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)

    return paths[0] if len(paths) == 1 else paths


####################################
##########  Sector plots ###########
####################################
def sector(r, theta, rlim=(0.0, 1.0), thetalim=(0.0, 360.0), clim=None, rotate=0.0, rlabel="", thetalabel="", clabel=None, label=None, rstep=None,
           thetastep=15.0, rticks='auto', thetaticks='auto', cbar_invert=False, fig=None, plot_kw={}, **kwargs):
    """ Sector Plot function

    Plots a sector plot (a.k.a "pizza plot") based on data with one radial axis and an angular axis

    Parameters
    ----------
    r : array-like or list
        Radial axis data.
    theta : array-like or list
        Angular axis data (degrees).
    rlim : tuple-like, optional
        The lower and upper limits for the radial axis (degrees).
    thetalim : tuple-like, optional
        The lower and upper limits for the angular axis (degrees).
    clim : tuple-like, optional
        Defines the limits of the colour-axis, it must contain two elements (lower and higer limits).
        Functions equivalently to the `vmin, vmax` arguments used by `colors.Normalize`. If both are
        given, `clim` takes priority.
    rotate : float, optional
        By how many degrees (clockwise) to rotate the entire plot (valid values in [-180, 180]).
    rlabel : str, optional
        Sets the label of the r-axis.
    thetalabel : str, optional
        Sets the label of the theta-axis.
    clabel : str, optional
        Setting `clabel` triggers the generation of a colourbar with axis label given by its value.
    label : str, optional
        Sets the label for the scatter plot.
    rstep : float, optional
        Sets the step size of r ticks.
    thetastep : float, optional, default: 15.0
        Sets the step size of theta ticks (degrees).
    rticks : 'auto', or ticker
        * Not implement *
    thetaticks : 'auto', or ticker
        * Not implement *
    cbar_invert : bool, optional
        If True, inverts the direction of the colour bar (not the colour map).
    fig : pyplot.Figure, optional
        Use the given figure to make the plot, defaults to the current figure.
    plot_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib scatter function.
        Parameters will be overwritten if also given implicitly in **kwargs.
    **kwargs : Collection properties, optional
        kwargs are used to specify matplotlib specific properties such as cmap, marker, norm, etc.
        A list of available `Collection` properties can be found here:
        https://matplotlib.org/3.1.0/api/collections_api.html#matplotlib.collections.Collection

    Returns
    -------
    ax : The pyplot.Axes object created for the sector plot.
    """

    from matplotlib.transforms import Affine2D
    from matplotlib.projections.polar import PolarAxes
    from matplotlib.pyplot import gcf, colorbar, legend

    from mpl_toolkits.axisartist import floating_axes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator, DictFormatter)
    import mpl_toolkits.axisartist.angle_helper as angle_helper

    from numpy import array, linspace, arange, shape, sqrt, floor, round, degrees, radians, pi

    if (fig is None):
        fig = gcf()

    # rotate a bit for better orientation
    trans_rotate = Affine2D().translate(0.0, 0)

    # scale degree to radians
    trans_scale = Affine2D().scale(pi / 180.0, 1.0)
    trans = trans_rotate + trans_scale + PolarAxes.PolarTransform()

    # Get theta ticks
    thetaticks = arange(*radians(array(thetalim) - rotate), step=radians(thetastep))
    theta_gridloc = FixedLocator(thetaticks[thetaticks / (2 * pi) < 1])
    theta_tickfmtr = DictFormatter(dict(zip(thetaticks, [f"{(round(degrees(tck)+rotate)):g}" for tck in thetaticks])))

    # tick_fmtr=DictFormatter(dict(angle_ticks))
    # tick_fmtr=angle_helper.Formatter()

    if (rstep is None):
        rstep = 0.5

    r_gridloc = FixedLocator(arange(rlim[0], rlim[1], step=rstep))

    grid = floating_axes.GridHelperCurveLinear(
        PolarAxes.PolarTransform(),
        extremes=(*radians(array(thetalim) - rotate), *rlim),
        grid_locator1=theta_gridloc,
        grid_locator2=r_gridloc,
        tick_formatter1=theta_tickfmtr,
        tick_formatter2=None,
    )

    ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid)
    fig.add_subplot(ax)

    # tick references
    thetadir_ref = ['top', 'right', 'bottom', 'left']
    rdir_ref = ['bottom', 'left', 'top', 'right']

    # adjust axes directions
    ax.axis["left"].set_axis_direction('bottom')  # Radius axis (displayed)
    ax.axis["right"].set_axis_direction('top')  # Radius axis (hidden)
    ax.axis["top"].set_axis_direction('bottom')  # Theta axis (outer)
    ax.axis["bottom"].set_axis_direction('top')  # Theta axis (inner)

    # Top theta axis
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction(thetadir_ref[(int(rotate) // 90) % 4])
    ax.axis["top"].label.set_axis_direction(thetadir_ref[(int(rotate) // 90) % 4])

    # Bottom theta axis
    ax.axis["bottom"].set_visible(False if rlim[0] < (rlim[1] - rlim[0]) / 3 else True)
    ax.axis["bottom"].major_ticklabels.set_axis_direction(thetadir_ref[(int(rotate) // 90 + 2) % 4])

    # Visible radius axis
    ax.axis["left"].major_ticklabels.set_axis_direction(rdir_ref[(int(rotate) // 90) % 4])
    ax.axis["left"].label.set_axis_direction(rdir_ref[(int(rotate) // 90) % 4])

    # Labels
    ax.axis["left"].label.set_text(rlabel)
    ax.axis["top"].label.set_text(thetalabel)

    # create a parasite axes whose transData in RA, cz
    sector_ax = ax.get_aux_axes(trans)

    # This has a side effect that the patch is drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to prevent this.
    sector_ax.patch = ax.patch
    sector_ax.patch.zorder = 0.9

    L = shape(theta)[0] if len(shape(theta)) > 1 else 1
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)

    # Insert clim as vmin, vmax into **kwargs dictionary, if given.
    if (clim is not None):
        try:
            _ = (e for e in clim)
            if (len(clim) == 2):
                plot_par['vmin'] = clim[0]
                plot_par['vmax'] = clim[1]
            else:
                raise TypeError("`clim` must be of iterable type and have two values only.")
        except (TypeError):
            raise TypeError("`clim` must be of iterable type and have two values only.")

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    # plot_par=dict_splicer(plot_par,L,[1]*L)

    if (L == 1):
        sctr = sector_ax.scatter(theta - rotate, r, label=label, **plot_par)
    else:
        for ii in range(L):
            sctr = sector_ax.scatter(theta[ii] - rotate, r[ii], label=label[ii], **plot_par[ii])

    if clabel is not None:
        cbar = colorbar(sctr)
        cbar.set_label(clabel)
        if cbar_invert:
            cbar.ax.invert_yaxis()

    return sector_ax


####################################
########  Statistics bands  ########
####################################
def statband(x, y, bin_type=None, bins=None, stat_mid='mean', stat_low='std', stat_high='std', from_mid=None, line=False, xlim=None, ylim=None,
             xinvert=False, yinvert=False, xlog=False, ylog=None, nmin=0, title=None, xlabel=None, ylabel=None,
             label=None, lab_loc=0, ax=None, grid=None, line_kw={}, band_kw={}, **kwargs):
    """Statistics line and band plotting function.

    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like.
    y : array-like or list
        If list it is assumed that each elemement is array-like.
    bin_type : {'number','width','edges','equal'}, optional
        Defines how is understood the value given in bins: 'number' for the desired number of bins, 'width' for the width
        of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements (or as close
        as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if float and 'edges'
        if ndarray.
    bins : int, float, array-like or list, optional
        Gives the values for the bins, according to bin_type.
    stat_mid : str, int, float or function, optional
        Defines how to calculate the midpoint of the statistics band. When passing a string it must be either one of the options
        for scipy.stats.binned_statistic(), i.e. 'mean', 'std', 'median', 'count', 'sum', 'min', 'max' or a user-defined function.
        If given as an integer or float, the number represents the value for the percentile to calculate in each bin.
        A function can be given which takes (only) a 1D array of values and returns a numerical statistic.
    stat_low / stat_high : str, int, float or function, optional
        Defines how to calculate the lower/upper limits for the statistic band. Can be given as one of the recognised strings above or as
        a string combining 'std' with a number, i.e. '[n]std', where [n] is the number of standard deviations away from the line of `stat_mid`.
        Can also be given as a number (integer or float) or function as described for stat_mid.
    from_mid : boolean, optional
        If True, the lower/upper bounds of the band are determined as the separation from the stat_mid line: i.e. stat_mid +/- stat_[low/high],
        otherwise, they are set to the values returned by stat_[low/high]. Defaults to True if stat_[low/high] are standard deviations.
    line : boolean, optional
        If True, draw a line that follows the statistic defined in line_stat.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If True, inverts the x-axis.
    yinvert : bool or list, optional
        If True, inverts the y-axis.
    xlog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    nmin : int, optional (default: 0)
        The minimum number of points required in a bin in order to be plotted.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the label for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, antialiasing, etc.
        A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    None
    """

    from splotch.base_func import axes_handler, bin_axis, plot_finalizer

    import numpy as np
    from numbers import Number
    import scipy.stats as stats
    from numpy import percentile
    from functools import partial
    import matplotlib.colors as clr
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import gca
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if ylog is None:
        from splotch.defaults import Params
        ylog = Params.hist1D_yaxis_log
    if bins is None:
        bins = int((len(x))**0.4)
    if 'linewidth' not in band_kw.keys():
        band_kw['linewidth'] = 0
    if 'alpha' not in band_kw.keys():
        band_kw['alpha'] = 0.4

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # band_par={**plot_kw, **kwargs} # For Python > 3.5
    band_kw.update(kwargs)

    # Check stat_low/stat_high arguments
    band_stat = np.array([None, None])
    band_multi = np.ones(2)
    for i, stat in enumerate([stat_low, stat_high]):  # loop over low/high statistic
        if isinstance(stat, Number):  # stat given as a percentile number
            band_stat[i] = partial(percentile, q=stat)  # Set as the percentile function with kwargs fixed.
        elif callable(stat):  # stat given as a function
            band_stat[i] = stat
        elif isinstance(stat, str) and 'std' in stat:  # stat given as a string with 'std'
            band_stat[i] = 'std'
            multi = list(filter(None, stat.split('std')))  # The multiplier for std, if any.
            if len(multi) == 0:  # No multiplier given
                band_multi[i] = 1.0
            elif len(multi) == 1:  # Multiplier was given
                band_multi[i] = multi[0]
            else:
                raise ValueError(f"Statistic '{stat}' not valid. Should be given as '[n]std', where [n] is a number.")
        else:
            raise ValueError(f"Statistic of type '{type(stat)}' was not recognised. Must be either a Number, function or string in the format '[n]std'.")

    # Check stat_mid argument
    if isinstance(stat_mid, Number):
        stat_mid = partial(percentile, q=stat_mid)

    # Assign 'from_mid' if not explicitly set
    if from_mid is None:
        if (band_stat[0] in ['std', np.std, np.nanstd] and (band_stat[1] in ['std', np.std, np.nanstd])):  # Band stats a type of standard deviation
            from_mid = True
        else:
            from_mid = False

    temp_x, bins_hist, bins_plot = bin_axis(x, bin_type, bins, log=xlog)
    temp_y = stats.binned_statistic(temp_x, y, statistic=stat_mid, bins=bins_hist)[0]
    if band_stat[0] == band_stat[1]:
        band_low, band_high = [stats.binned_statistic(temp_x, y, statistic=band_stat[0], bins=bins_hist)[0]] * 2
    else:
        band_low = stats.binned_statistic(temp_x, y, statistic=band_stat[0], bins=bins_hist)[0]
        band_high = stats.binned_statistic(temp_x, y, statistic=band_stat[1], bins=bins_hist)[0]

    if from_mid is True:  # Band intervals should be taken as the difference from the mid line
        band_low = temp_y - band_multi[0] * band_low
        band_high = temp_y + band_multi[1] * band_high

    if ylog:
        temp_y = np.where(temp_y == 0, np.nan, temp_y)
    
    counts = np.histogram(temp_x, bins=bins_hist)[0]
    x = stats.binned_statistic(temp_x, x, statistic=stat_mid, bins=bins_hist)[0]
    y = temp_y
    
    x = x[counts>nmin]
    y = y[counts>nmin]
    band_low = band_low[counts>nmin]
    band_high = band_high[counts>nmin]
    
    plt.fill_between(x, band_low, band_high, label=label, **band_kw)

    if line:
        plt.plot(x, y, **line_kw)
    if label is not None:
        plt.legend(loc=lab_loc)

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)


####################################
########  Statistics bars  #########
####################################
def statbar(x, y, bin_type=None, bins=None, stat_cen='mean', bar_x=True, stat_y='std', line=False, xlim=None, ylim=None,
            xinvert=False, yinvert=False, xlog=False, ylog=None, title=None, xlabel=None, ylabel=None,#nmin=0,
            label=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):
    """Statistics line and bar plotting function.

    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like.
    y : array-like or list
        If list it is assumed that each elemement is array-like.
    bin_type : {'number','width','edges','equal'}, optional
        Defines how is understood the value given in bins: 'number' for the desired number of bins, 'width' for the width
        of the bins, 'edges' for the edges of bins, and 'equal' for making bins with equal number of elements (or as close
        as possible). If not given it is inferred from the data type of bins: 'number' if int, 'width' if float and 'edges'
        if ndarray.
    bins : int, float, array-like or list, optional
        Gives the values for the bins, according to bin_type.
    stat_cen : str, int, float, function, or 2-element array-like of any of the previous type, optional
        Defines how to calculate the position of centre of each errorbar. When passing an integer or float is
        interpreted as being the percentile for the limit. When passing a function it must have the input and
        ouput characteristics required by scipy.stats.binned_statistic().
    bar_x : bool, optional
        If False turns off the display of the bin widths with bars.
    stat_y : str, int, float, function, or 2-element array-like of any of the previous type, optional
        Defines how to calculate the y error bars. When passing a string it must be either one of the options
        for scipy.stats.binned_statistic(), or a string that combines 'std' with a number (e.g., '2.2std'), where to number is
        interpreted as the number of standard deviations that the limit must cover. When passing an integer or float is
        interpreted as being the percentile for the limit. When passing a function it must have the input and ouput
        characteristics required by scipy.stats.binned_statistic().
    line : boolean, optional
        If True, draw a line that follows the statistic defined in line_stat.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If True, inverts the x-axis.
    yinvert : bool or list, optional
        If True, inverts the y-axis.
    xlog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True, the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the label for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth, antialiasing, etc.
        A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    None
    """

    from splotch.base_func import axes_handler, bin_axis, plot_finalizer

    import numpy as np
    from numbers import Number
    import scipy.stats as stats
    from numpy import percentile
    from functools import partial
    import matplotlib.colors as clr
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import gca, errorbar, rcParams
    
    # Set the current axis
    if ax is not None:
        old_axes=axes_handler(ax)
    else:
        ax=gca()
        old_axes=ax
    
    if ylog is None:
        from splotch.defaults import Params
        ylog = Params.hist1D_yaxis_log
    if bins is None:
        bins = int((len(x))**0.4)

    if not isinstance(stat_y, str):
        try:
            iter(stat_y)
        except TypeError:
            stat_y = [stat_y, stat_y]
    else:
        stat_y = [stat_y, stat_y]

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par = {**plot_kw, **kwargs}  # For Python > 3.5
    plot_par = plot_kw.copy()
    plot_par.update(kwargs)
    if 'linewidth' not in plot_par.keys():
        plot_par['linewidth'] = 0
    if 'elinewidth' not in plot_par.keys():
        plot_par['elinewidth'] = rcParams['lines.linewidth']

    bar_multi = np.ones(2)
    for i in range(len(stat_y)):
        if isinstance(stat_y[i], Number):
            stat_y[i] = partial(percentile, q=stat_y[i])
        elif 'std' in stat_y[i] and len(stat_y[i].replace('std', '')) > 0:
            bar_multi[i] = float(stat_y[i].replace('std', ''))
            stat_y[i] = 'std'

    if isinstance(stat_cen, Number):
        stat_cen = partial(percentile, q=stat_cen)

    temp_x, bins_hist, bins_plot = bin_axis(x, bin_type, bins, log=xlog)
    temp_y = stats.binned_statistic(temp_x, y, statistic=stat_cen, bins=bins_hist)[0]
    if stat_y[0] == stat_y[1]:
        bar_low, bar_high = [stats.binned_statistic(temp_x, y, statistic=stat_y[0], bins=bins_hist)[0]] * 2
    else:
        bar_low = stats.binned_statistic(temp_x, y, statistic=stat_y[0], bins=bins_hist)[0]
        bar_high = stats.binned_statistic(temp_x, y, statistic=stat_y[1], bins=bins_hist)[0]
    if stat_y[0] == 'std':
        bar_low = bar_multi[0] * bar_low
    else:
        bar_low = temp_y - bar_low
    if stat_y[0] == 'std':
        bar_high = bar_multi[0] * bar_high
    else:
        bar_high -= temp_y
    if ylog:
        temp_y = np.where(temp_y == 0, np.nan, temp_y)

    x = stats.binned_statistic(temp_x, x, statistic=stat_cen, bins=bins_hist)[0]
    y = temp_y
    if bar_x:
        bar_x = [x - bins_plot[:-1], bins_plot[1:] - x]

    if bar_x:
        errorbar(x, y, xerr=bar_x, yerr=[bar_low, bar_high], **plot_par)
    else:
        errorbar(x, y, yerr=[bar_low, bar_high], **plot_par)

    if label is not None:
        plt.legend(loc=lab_loc)

    plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, grid)
    if old_axes is not ax:
        old_axes = axes_handler(old_axes)
