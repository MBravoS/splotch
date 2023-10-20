########################################################################
############## Definition of all wrappers for 1D plotting ##############
########################################################################
from collections.abc import Iterable
from warnings import warn
from numpy import array, diff, dtype, empty, full_like, histogram, inf, linspace, log10, logspace, meshgrid, nan, nancumsum, arange
from numpy import nanmax, nanmean, nanmin, nanstd, nansum, ndarray, ones, piecewise, prod, shape, squeeze, where
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from matplotlib.legend import Legend
from matplotlib.legend_handler import HandlerLine2D, HandlerPathCollection, HandlerTuple
from matplotlib.pyplot import bar, barh, fill_between, fill_betweenx, gcf, gca, legend, plot as plt_plot, sca, stairs  # step
from matplotlib.transforms import Bbox
from matplotlib.colors import to_rgba_array
from scipy.stats import binned_statistic
from sympy import Expr, latex, sympify
from sympy.utilities.lambdify import lambdify
import inspect

from .base_func import axes_handler, grid_handler, lims_handler, bin_axis, dict_splicer, is_numeric, is_number, is_listlike, plot_finalizer, _plot_finalizer, simpler_dict_splicer
from .defaults import Params


####################################
# Generalized lines
####################################
def axline(x=None, y=None, a=None, b=None,
           xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False, title=None,
           xlabel=None, ylabel=None, label=None, grid=None, ax=None, plot_kw={}, **kwargs):
    """Generalised axis lines.

    This function aims to generalise the usage of axis lines calls (axvline/axhline) together and
    to allow lines to be specified by a slope/intercept according to the function y=a*x + b.

    Parameters
    ----------
    x : int or list, optional
        x position(s) in data coordinates for a vertical line(s).
    y : int or list, optional
        y position(s) in data coordinates for a horizontal line(s).
    a : int or list, optional
        Slope(s) of diagonal axis line(s), defaults to 1 if not specified when b is given.
    b : int or list, optional
        Intercept points(s) of diagonal axis line(s), defaults to 0 if not specified when a is given.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If true inverts the x-axis.
    yinvert : bool or list, optional
        If true inverts the y-axis.
    xlog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the legend for the plot.
    grid : boolean, dict or None
        The grid behaviour, acts according to:
         * If boolean, the grid is set on (True) or off (False).
         * If dict, allows specific `.Line2D` properties of the grid to be set.
         * If None, use default grid parameters if this is the initial plotting call, otherwise do nothing.
    ax : pyplot.Axes or list-like, optional
        Specifies the axes on which to plot the lines, defaults to the current axes.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
        antialiasing, etc. A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    lines
        A list of Line2D objects representing the plotted lines. If more than one subplot and/or line is plotted,
        the resulting list will be multi-dimensional with the dimensions corresponding to the axes/lines.

    """

    # Get the relevant axis
    if ax is not None:
        if is_listlike(ax):
            if len(shape(ax)) > 1:  # If ax array is multi-dimensional, flatten it
                ax = array(ax).flatten()
        else:
            ax = [ax]  # Axis must be a list-like object
    else:
        ax = [gca()]
    
    old_ax = axes_handler(ax[0])  # sets the current axis and returns old axis

    # Validate input parameters
    if not (any([is_numeric(var) for var in [x, y, a, b]])):  # If nothing has been specified
        raise TypeError("axline() missing one of optional arguments: 'x', 'y', 'a' or 'b'")

    # If parameter given and it is not a list-like object, convert to list
    if x and not is_listlike(x):
        x = [x]

    if y and not is_listlike(y):
        y = [y]

    if a and not is_listlike(a):
        a = [a]

    if b and not is_listlike(b):
        b = [b]

    if (x is not None and y is not None):  # Check whether both x and y were specified
        raise ValueError("'x' and 'y' cannot be both specified")

    if (x is not None):  # Check conditions if x specified
        if (any([a, b])):  # Should not specify a or b, if x given.
            raise ValueError("'{0}' cannot be specified if x specified".format('a' if a else 'b'))
        L = len(x)

    if (y is not None):  # Check conditions if y specified
        if (any([a, b])):  # Should not specify a or b, if y given.
            raise ValueError("'{0}' cannot be specified if y specified".format('a' if a else 'b'))
        L = len(y)

    if (a is not None):
        if (b is None):  # If no intercept specified
            b = [0] * len(a)  # set b to 0 for all a
        else:
            if (len(b) == 1):
                b = [b[0]] * len(a)
            elif (len(b) != len(a)):
                if (len(a) == 1):
                    a = [a[0]] * len(b)
                else:
                    raise ValueError(f"Length of 'a' ({len(a)}) and length of 'b' ({len(b)}) must be equal or otherwise 1")
        L = len(a)
    elif (b is not None):
        if (a is None):  # If no slope specified
            a = [1] * len(b)  # set a to 1 for all b
        L = len(b)

    if not isinstance(label, list):
        label = [label] * L
    elif (len(label) != L):
        raise ValueError("Length of label list ({0}) must match the number of lines given ({1}).".format(len(label), L))

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    plot_par = {**plot_kw, **kwargs}

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    lines = [[]] * len(ax)  # Initialise list which contains each Line2D object
    for jj, axis in enumerate(ax):  # Loop over all axes
        gridpar = grid_handler(grid, axis)

        ax_xlim = axis.get_xlim() if lims_handler(xlim, axis) is None else xlim
        ax_ylim = axis.get_ylim() if lims_handler(ylim, axis) is None else ylim

        if (x is not None):
            for ii, xx in enumerate(x):
                lines[jj].append(axis.axvline(x=xx, **plot_par[ii], label=label[ii]))
        if (y is not None):
            for ii, yy in enumerate(y):
                lines[jj].append(axis.axhline(y=yy, **plot_par[ii], label=label[ii]))
        if (a is not None):
            for ii, (aa, bb) in enumerate(zip(a, b)):  # Loop over all lines
                lines[jj].append(axis.axline(xy1=(0, bb), slope=aa, label=label[ii], **plot_par[ii]))

        _plot_finalizer(xlog, ylog, ax_xlim, ax_ylim, title, xlabel, ylabel, xinvert, yinvert, gridpar, axis)

    if old_ax is not None:  # Reset the previously set axis
        sca(old_ax)

    return squeeze(lines).tolist()  # Reduce the dimensionality of the lines, if needed


####################################
# Broken axis plot
####################################
def brokenplot(x, y=None, xbreak=None, ybreak=None, xlim=None, ylim=None, sep=0.05, overflow=0.05,
               xinvert=False, yinvert=False, xlog=False, ylog=False, title=None,
               xlabel=None, ylabel=None, label=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):
    """Broken Axis Plot Function
    
    Creates a standard plot call with an axis break at `xbreak` or `ybreak` for vertical or
    horizontal breaks.
    
    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like. If y is not given, the given values
        pass to y and a numpy array is generated with numpy.arange() for the x values.
    y : array-like or list, optional
        If list it is assumed that each elemement is array-like.
    xbreak/ybreak : float or tuple-like, required
        The location(s) of the vertical or horizontal breaks is controlled by xbreak or ybreak, respectively.
        The value can be a single location or a tuple defining the (start, stop) points of the break.
        Only one coordinate can be broken in a given plot.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    sep : float, optional (default: 0.05)
        The separation size of the axis break, given as a fraction of the axis dimensions.
    overflow : float, optional (default: 0.05)
        The fractional overflow of the axes beyond the break values. Often good to keep this non-zero as
        breaking the axis exactly at the break point leaves the final gridlines and ticks out of view.
    xinvert : bool or list, optional
        If true inverts the x-axis.
    yinvert : bool or list, optional
        If true inverts the y-axis.
    xlog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the legend for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, optional
        If not given defaults to the value defined in splotch.Params.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
        antialiasing, etc. A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
    
    Returns
    -------
    lines : list of Line2D
        A list of Line2D objects (paired as tuples) representing the plotted data.
        The lines are given as pairs to correspond to the separate lines either side of the x/ybreak.
    ax2 : pyplot.Axes, optional
        The additional axis object created.
    """
    
    # Get the relevant axis
    ax1 = gca() if ax is None else ax
    old_ax = axes_handler(ax1)

    gridpar = grid_handler(grid, ax1)

    if not is_listlike(x) or len(shape(x)) == 1:
        x = [x]

    L = len(x)

    if y is None:
        y = x
        x = [arange(len(x[i])) for i in range(L)]
    else:
        if not is_listlike(y) or len(shape(y)) == 1:
            y = [y]

    if not is_listlike(label) or len(shape(label)) == 1:
        label = [label] * L
        
    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    plot_par = {**plot_kw, **kwargs}
    
    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)
    
    # Get the original axis position
    pos0 = ax.get_position(original=True)
    width0, height0 = pos0.x1 - pos0.x0, pos0.y1 - pos0.y0
    
    # Perform first plot call to original axis
    lines = []
    for i in range(L):
        lines += ax1.plot(x[i], y[i], label=label[i], **plot_par[i])
    
    # Get the axis limits if not already specified
    xlims = ax1.get_xlim() if xlim in [None, 'auto'] else xlim
    ylims = ax1.get_ylim() if ylim in [None, 'auto'] else ylim
    
    # Validate xbreak/ybreak
    if (xbreak is None and ybreak is None):
        raise ValueError("Require either xbreak/ybreak to be specified.")

    if (ybreak is not None and xbreak is not None):
        raise ValueError("Cannot specify both xbreak and ybreak.")

    breaks = xbreak if xbreak else ybreak
    if not is_listlike(breaks):
        breaks = (breaks,) * 2
    else:
        if len(breaks) == 1:
            breaks = (breaks[0],) * 2
        elif len(breaks) == 2:
            if breaks[0] > breaks[1]:
                raise ValueError("First value in xbreak/ybreak must not be greater than the second value.")
        else:
            raise ValueError("xbreak/ybreak must be a single value or a list-like object with two elements.")

    if xbreak and (breaks[0] < xlims[0] or breaks[1] > xlims[1]):
        raise ValueError(f"xbreak given ({xbreak}) must be within the xlim ({xlims}).")
    if ybreak and (breaks[0] < ylims[0] or breaks[1] > ylims[1]):
        raise ValueError(f"ybreak given ({ybreak}) must be within the ylim ({ylims}).")

    # Define the positions of the two separated axes
    pos1 = Bbox(list(pos0.get_points()))
    if xbreak:
        pos1.x1 = pos1.x0 + (pos1.x1 - pos1.x0) * (sum(breaks) / 2 - xlims[0]) / (xlims[1] - xlims[0]) - sep * (pos1.x1 - pos1.x0) / 2
    if ybreak:
        pos1.y1 = pos1.y0 + (pos1.y1 - pos1.y0) * (sum(breaks) / 2 - ylims[0]) / (ylims[1] - ylims[0]) - sep * (pos1.y1 - pos1.y0) / 2

    pos2 = Bbox(list(pos0.get_points()))
    if xbreak:
        pos2.x0 = pos2.x0 + (pos2.x1 - pos2.x0) * (sum(breaks) / 2 - xlims[0]) / (xlims[1] - xlims[0]) + sep * (pos2.x1 - pos2.x0) / 2
    if ybreak:
        pos2.y0 = pos2.y0 + (pos2.y1 - pos2.y0) * (sum(breaks) / 2 - ylims[0]) / (ylims[1] - ylims[0]) + sep * (pos2.y1 - pos2.y0) / 2
    
    ax1.set_position(pos1)  # Resize the first axis
    ax2 = ax1.figure.add_axes(pos2)  # Add and duplicate the plotting in the second axis
    
    # Perform first plot call to original axis
    lines = []
    for i in range(L):
        lines += ax2.plot(x[i], y[i], label=label[i], **plot_par[i])
    
    width1, height1 = pos1.x1 - pos1.x0, pos1.y1 - pos1.y0
    width2, height2 = pos2.x1 - pos2.x0, pos2.y1 - pos2.y0

    # Adjust x/y limits for both axes
    xlim1 = (xlims[0], breaks[0] + overflow * (width1 + sep) / width0) if xbreak else xlims
    xlim2 = (breaks[1] - overflow * (width2 + sep) / width0, xlims[1]) if xbreak else xlims
    
    ylim1 = (ylims[0], breaks[0]) if ybreak else ylims
    ylim2 = (breaks[1], ylims[1]) if ybreak else ylims
    
    if xbreak:
        dx1, dy1 = 0.01 * width0 / (width0 - width2 - sep / 2), height1 * 0.025
        dx2, dy2 = 0.01 * width0 / (width0 - width1 - sep / 2), height2 * 0.025
        
        dash_kw = dict(transform=ax1.transAxes, color='black', linestyle='-', marker='', clip_on=False)
        temp1 = (breaks[0] - xlim1[0]) / (breaks[0] - xlim1[0] + overflow * (width1 + sep) / width0)
        ax1.plot((temp1 - dx1, temp1 + dx1), (0 - dy1, 0 + dy1), **dash_kw)  # bottom-left
        ax1.plot((temp1 - dx1, temp1 + dx1), (1 - dy1, 1 + dy1), **dash_kw)  # top-left
        
        dash_kw.update(transform=ax2.transAxes)  # switch to the left axes
        temp2 = 1 - (xlim2[1] - breaks[1]) / (xlim2[1] - breaks[1] + overflow * (width2 + sep) / width0)
        ax2.plot((temp2 - dx2, temp2 + dx2), (0 - dy2, 0 + dy2), **dash_kw)  # bottom-right
        ax2.plot((temp2 - dx2, temp2 + dx2), (1 - dy2, 1 + dy2), **dash_kw)  # top-right

        # Fix axes labels and splines
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(labelright=False, which='both')
        ax1.yaxis.tick_left()

        ax2.spines['left'].set_visible(False)
        ax2.tick_params(labelleft=False, which='both')
        ax2.yaxis.tick_right()

        if (ax1.get_xticks()[-1] == ax2.get_xticks()[0]):
            if (breaks[0] >= (xlims[0] + xlims[1]) * 0.5):
                ax1.set_xticks(ax1.get_xticks()[:-1])  # Remove duplicate tick on left side
            else:
                ax2.set_xticks(ax2.get_xticks()[1:])  # Remove duplicate tick on right side

        ax1.axes.spines[['top', 'bottom']].set_bounds((xlim1[0], breaks[0]))
        ax2.axes.spines[['top', 'bottom']].set_bounds((breaks[1], xlim2[-1]))
        
        ax1.xaxis.set_label_coords(0.5 * (width0 / width1), -0.125, transform=ax1.axes.transAxes)
                
    if ybreak:
        dx1, dy1 = width1 * 0.0125, 0.02 * height0 / (height0 - height2 - sep / 2)
        dx2, dy2 = width2 * 0.0125, 0.02 * height0 / (height0 - height1 - sep / 2)
        
        dash_kw = dict(transform=ax1.transAxes, color='black', linestyle='-', marker='', clip_on=False)
        temp1 = (breaks[0] - ylim1[0]) / (breaks[0] - ylim1[0] + overflow * (height1 + sep) / height0)
        ax1.plot((0 - dx1, 0 + dx1), (temp1 - dy1, temp1 + dy1), **dash_kw)  # bottom-left
        ax1.plot((1 - dx1, 1 + dx1), (temp1 - dy1, temp1 + dy1), **dash_kw)  # bottom-right
        
        dash_kw.update(transform=ax2.transAxes)  # switch to the left axes
        temp2 = 1 - (ylim2[1] - breaks[1]) / (ylim2[1] - breaks[1] + overflow * (height2 + sep) / height0)
        ax2.plot((0 - dx2, 0 + dx2), (temp2 - dy2, temp2 + dy2), **dash_kw)  # top-left
        ax2.plot((1 - dx2, 1 + dx2), (temp2 - dy2, temp2 + dy2), **dash_kw)  # top-right

        # Fix axes labels and splines
        ax1.spines['top'].set_visible(False)
        ax1.tick_params(labeltop=False, which='both')
        ax1.xaxis.tick_bottom()

        ax2.spines['bottom'].set_visible(False)
        ax2.tick_params(labelbottom=False, which='both')
        ax2.xaxis.tick_top()

        if (ax1.get_yticks()[-1] == ax2.get_yticks()[0]):
            if (breaks[0] >= (ylims[0] + ylims[1]) * 0.5):
                ax1.set_yticks(ax1.get_yticks()[:-1])  # Remove duplicate tick on top
            else:
                ax2.set_yticks(ax2.get_yticks()[1:])  # Remove duplicate tick on bottom
                
        ax1.axes.spines[['left', 'right']].set_bounds((ylim1[0], breaks[0]))
        ax2.axes.spines[['left', 'right']].set_bounds((breaks[1], ylim2[-1]))

        ax1.yaxis.set_label_coords(-0.125, 0.5 * (height0 / height1), transform=ax1.axes.transAxes)
    
    if any(label):
        ax2.legend(loc=lab_loc)
    
    _plot_finalizer(xlog, ylog, xlim1, ylim1, title, xlabel, ylabel, xinvert, yinvert, gridpar, ax=ax1)
    _plot_finalizer(xlog, ylog, xlim2, ylim2, None, None, None, xinvert, yinvert, gridpar, ax=ax2)
    
    # Share axes
    if xbreak:
        ax1.get_shared_y_axes().join(ax, ax2)
    if ybreak:
        ax1.get_shared_x_axes().join(ax, ax2)
    
    if old_ax is not None:
        sca(old_ax)

    return (lines[0] if len(lines) == 1 else lines, ax2)


########################################
# Curves from mathematical expressions #
########################################
def curve(expr, var=None, subs={}, orientation='horizontal', permute=False, bounds=None, num=101,
          xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False, grid=None,
          title=None, xlabel=None, ylabel=None, label=None, uselatex=True, ax=None, plot_kw={}, **kwargs):
    """Plot Mathematical Expressions
    
    Plot the curve(s) corresponding to mathematical expressions over a given range across the independent variable.
    Expressions can be given with multiple variables, whereby one is taken to be the independent variable and all
    others are substitution variables, which can take a variable number of values producing the corresponding number
    of curves.
    
    Parameters
    ----------
    expr : str, sympy.Expr, callable() or list-like
        An expression parsed either as a string, sympy expression or callable (function or lambda)
        which will be evaluated by the function in the range of `bounds`.
    var : str or sympy symbol, default: 'x'
        The independent variable on which to evaluate the expression (i.e. the variable on the x-axis).
        This defaults to the first non-numeric element of the expression or otherwise simply assumes
        this to be 'x' (if `orientation='horizotnal'`, else 'y').
    subs : dict, optional
        If `expr` contains more symbols than the independent variable `var`, this dictionary will
        substitute numerical values for all additonal symbols given. `subs` is required if additional
        symbols are specified in the expression.
    permute : bool, optional (default: False)
        If `subs` contains more than one substitution variable with multiple values, setting `permute=True`
        will generate curves for each permutation of values given.
        For example, `splotch.curve('a*x^2 + b',subs = {'a'=[1,2],'b'=[3,4,5]})` will create six curves.
    orientation : str, optional (default: 'horizontal')
        The orientation of the independent axis, i.e. whether the independent variable is defined along
        the x-axis ('horizontal') or the y-axis ('vertical') of the plot.
        splotch.curve('a*x') is notationally the same as splotch.curve('1/a*y',var='y',orientation='vertical').
    bounds : list-like, optional
        The range over which the function will be plotted. If not given, these default to
        the current bounds of the plot.
    num : int, optional (default: 101)
        The number of values along the independent variable on which to evaulate `expr`.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If true inverts the x-axis.
    yinvert : bool or list, optional
        If true inverts the y-axis.
    xlog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    grid : boolean, dict or None
        The grid behaviour, acts according to:
         * If boolean, the grid is set on (True) or off (False).
         * If dict, allows specific `.Line2D` properties of the grid to be set.
         * If None, use default grid parameters if this is the initial plotting call, otherwise do nothing.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : bool, str or list-like, optional (default: True if `expr` is not callable)
        Sets the label(s) of the curve(s) for the legend. `label` can be one of:
            - `True`:
                Creates a label for every curve defined by `subs`. The label will list all values
                `subs` that produced the curve. If a parameter in `subs` has only one value
                (i.e. constant amongst all curves), it will not appear in the label.
            - `str`:
                If a single string is given, only one label will be shown and all curves will
                be shown as the label handle.
            - list-like (length must equal to number of curves):
                A label given to each curve that is produced.
            - `False` or `None`:
                No label will be assigned and a legend will not be shown.
    uselatex : bool, optional (default: True)
        If `label==True`, sets whether to use LaTeX when creating the labels for legend.
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D
        properties. It is recommended that kwargs be parsed implicitly through **kwargs
        for readability.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
        antialiasing, etc. A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    curves : list of (or single) pyplot.Line2D object(s)
        A list of Line2D objects created for each curved create by `subs`.
    expr : Sympy.Expr
        If expr was given as a string, this returns the sympy expression created from `sympy.sympify()`.
        Otherwise, simply returns the `expr` that was given.
    
    """
    
    # Get the relevant axis
    if ax is None: ax = gca()
    old_ax = axes_handler(ax)

    gridpar = grid_handler(grid, ax)
    xlim = lims_handler(xlim, ax)
    ylim = lims_handler(ylim, ax)

    # Assign bounds if none given
    if (bounds is None):
        if orientation == 'horizontal':
            bounds = xlim if xlim not in [None, 'auto'] else ax.get_xlim()
        else:
            bounds = ylim if ylim not in [None, 'auto'] else ax.get_ylim()
    
    # Parse expression
    isfunc = False
    if (isinstance(expr, str)):
        expr = sympify(expr)
    elif (callable(expr)):
        isfunc = True
    elif (isinstance(expr, Expr)):
        pass
    else:
        raise TypeError(f"`expr` must be of type `str`, sympy.Expr or callable, instead got {type(expr)}.")

    if not isfunc:  # expr is a Sympy expression
        symbols = expr.free_symbols  # Get all of the Symbols in the expression
        symbolkeys = [str(symb) for symb in symbols]  # convert these to strings, instead of sympy.Symbols
        
        if var is None:  # Assume independent variable is 'x', otherwise, assume the first symbol.
            if orientation == 'horizontal':
                var = 'x'  # if 'x' in symbolkeys or len(symbolkeys)==0 else None
            else:  # first test for 'y' as an independent variable, then default to x.
                if 'y' in symbolkeys:
                    var = 'y'
                else:
                    var = 'x'  # if 'x' in symbolkeys or len(symbolkeys)==0 else symbolkeys[0]
        
        # Validate the substitution variable names
        if subs is None:
            subs = dict()

        if (var in list(subs)):
            raise ValueError(f"Independent variable '{var}' should not be in subs")
            
        for key in list(subs):
            if (key not in symbolkeys):
                raise KeyError(f"Substitution variable '{key}' does not exist in 'expr'")
        
        # Check for missing substitution variables
        missing = set(symbolkeys) - set(subs) - set(var)
        if missing:
            raise TypeError(f"`expr` missing {len(missing)} required substitution variable{'s' if len(missing) > 1 else ''}: {list(missing)}")
    
    # The lengths of each substitute value list, len=1 if just a single value
    lens = [len(subs[key]) if (isinstance(subs[key], Iterable) and not isinstance(subs[key], str)) else 1 for key in list(subs)]
    if (permute is True):
        L = prod(lens)
        perms = array(meshgrid(*subs.values())).reshape(len(subs), -1)
        permsubs = {}
        for ii, key in enumerate(list(subs)):
            permsubs[key] = perms[ii]
        subsarr = simpler_dict_splicer(permsubs, L, [1] * L)
    else:
        L = max(lens) if len(lens) > 0 else 1
        subsarr = simpler_dict_splicer(subs, L, [1] * L)
    
    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    plot_par = {**plot_kw, **kwargs}
    
    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    # Set default label behaiviour
    if label is None:
        label = False if isfunc else True

    # Create the labels
    if label is False:
        labellist = [None] * L
    elif label is True:  # Auto-generate labels
        if isfunc:
            exprstr = expr.__name__
            signature = inspect.signature(expr)  # get the signature of the function/lambda
            pars = [str(par) for par in signature.parameters.values()]
            
            labellist = []
            for ii in range(L):  # Make a label for each of sub values
                labellist.append(f"{exprstr}(" + ", ".join(f"{par}={subsarr[ii][par]}" if par in list(subsarr[ii]) else f"{par}" for par in pars) + ")")

        elif subsarr == [{}]:
            labellist = [f"${latex(expr)}$" if uselatex else str(expr)]
        else:
            labellist = []
            exprstr = f"${latex(expr)}$" if uselatex else str(expr)
            for ii in range(L):  # Make a label for each of sub values
                if uselatex:
                    labellist.append(f"{exprstr} {{" + ", ".join([f"${key}$={subsarr[ii][key]}" for key in list(subsarr[ii])]) + "}")
                else:
                    labellist.append(f"{exprstr} {{" + ", ".join([f"{key}={subsarr[ii][key]}" for key in list(subsarr[ii])]) + "}")
    elif isinstance(label, str):  # A single string
        labellist = [label] * L
    else:
        if is_listlike(label):
            if (len(label) != L):
                raise TypeError(f"Number of labels ({len(label)}) does not match the number of curves ({L}).")
            else:
                labellist = label
        else:  # was not an iterable
            raise TypeError(f"`label` of type {type(label)} is not recognised.")

    # Create and plot the curves
    vararr = logspace(*log10(bounds), num=num) if xlog else linspace(*bounds, num=num)
    curves = [None] * L
    for ii in range(L):
        if (isfunc):
            curvearr = expr(vararr, **subsarr[ii])
        else:
            lamb = lambdify(var, expr.subs(subsarr[ii]), modules='numpy')  # returns a numpy-ready function

            if expr.subs(subsarr[ii]).is_constant():
                func = lambda x: full_like(x, lamb(x))
            else:
                func = lambda x: lamb(array(x))
            
            curvearr = func(vararr)
            
        curves[ii] = plt_plot(vararr if orientation == 'horizontal' else curvearr,
                              curvearr if orientation == 'horizontal' else vararr,
                              label=labellist[ii], **plot_par[ii])[0]
    
    _plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, gridpar, ax)
    
    if old_ax is not None:
        sca(old_ax)
    
    return(curves[0] if len(curves) == 1 else curves, expr)


##################################################
# Piecewise curves from mathematical expressions #
##################################################
def curve_piecewise(expr, var=None, subs={}, orientation='horizontal', bounds=None, intervals=[], permute=False, num=101,
                    xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False, grid=None,
                    title=None, xlabel=None, ylabel=None, label=None, uselatex=True, ax=None, plot_kw={}, **kwargs):
    """Plot Piecewise Mathematical Expressions
    
    Plot the curve(s) corresponding to mathematical expressions over a given range across the independent variable.
    Expressions can be given with multiple variables, whereby one is taken to be the independent variable and all
    others are substitution variables. This function can only accept one value per substitution variable.
    
    Parameters
    ----------
    expr : str, sympy.Expr, callable() or list-like
        An expression parsed either as a string, sympy expression or callable (function or lambda)
        which will be evaluated by the function in the range of `bounds`. A piece-wise function is defined
        by parsing a list of expressions. The piece-wise functionality must also be reflected in `bounds`.
    var : str or sympy symbol, required.
        The independent variable on which to evaluate the expression (i.e. the variable on the x-axis).
    subs : dict, optional
        If `expr` contains more symbols than the independent variable `var`, this dictionary will
        substitute numerical values for all additonal symbols given. `subs` is required if additional
        symbols are specified in the expression.
    orientation : str, optional (default: 'horizontal')
        The orientation of the independent axis, i.e. whether the independent variable is defined along
        the x-axis ('horizontal') or the y-axis ('vertical') of the plot.
        splotch.curve('a*x') is notationally the same as splotch.curve('1/a*y',var='y',orientation='vertical').
    bounds : list-like, optional
        The range over which the function will be plotted. If not given, these default to the current bounds
        of the axes being plotted onto, given by `ax`.
    intervals : list-like, required
        The points at which the curve partitions each of the N expressions given in `expr`. The values must lie
        lie within `bounds` as these are assumed to be the minimum and maximum points of the intervals.
    num : int, optional (default: 101)
        The number of values along the independent variable on which to evaulate `expr`.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If true inverts the x-axis.
    yinvert : bool or list, optional
        If true inverts the y-axis.
    xlog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    grid : boolean, dict or None
        The grid behaviour, acts according to:
         * If boolean, the grid is set on (True) or off (False).
         * If dict, allows specific `.Line2D` properties of the grid to be set.
         * If None, use default grid parameters if this is the initial plotting call, otherwise do nothing.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : bool, str or list-like, optional (default: True if `expr` is not callable)
        Sets the label(s) of the curve(s) for the legend. `label` can be one of:
            - `True`:
                Creates a label for every curve defined by `subs`. The label will list all values
                `subs` that produced the curve. If a parameter in `subs` has only one value
                (i.e. constant amongst all curves), it will not appear in the label.
            - `str`:
                If a single string is given, only one label will be shown and all curves will
                be shown as the label handle.
            - list-like (length must equal to number of curves):
                A label given to each curve that is produced.
            - `False` or `None`:
                No label will be assigned and a legend will not be shown.
    uselatex : bool, optional (default: True)
        If `label==True`, sets whether to use LaTeX when creating the labels for legend.
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are Line2D
        properties. It is recommended that kwargs be parsed implicitly through **kwargs
        for readability.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
        antialiasing, etc. A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    curves : list of (or single) pyplot.Line2D object(s)
        A list of Line2D objects created for each curved create by `subs`.
    expr : Sympy.Expr
        If expr was given as a string, this returns the sympy expression created from `sympy.sympify()`.
        Otherwise, simply returns the `expr` that was given.
    """

    # Get the relevant axis
    if ax is None: ax = gca()
    old_ax = axes_handler(ax)

    gridpar = grid_handler(grid, ax)
    xlim = lims_handler(xlim, ax)
    ylim = lims_handler(ylim, ax)
    
    # Assign bounds if none given
    if (bounds is None):
        if orientation == 'horizontal':
            bounds = xlim if xlim is not None else ax.get_xlim()
        else:
            bounds = ylim if ylim is not None else ax.get_ylim()

    # Check if iterable
    if is_listlike(expr):
        if isinstance(expr, dict):
            exprs = list(expr.values())
        else:
            exprs = list(expr)
    else:  # a string or something else
        exprs = [expr]

    # Parse expressions
    isfunc = [False] * len(exprs)
    
    for ii in range(len(exprs)):
        if (isinstance(exprs[ii], str)):
            exprs[ii] = sympify(exprs[ii])
        elif (callable(exprs[ii])):
            isfunc[ii] = True
        elif (isinstance(exprs[ii], Expr)):
            pass
        else:
            raise TypeError(f"Elements of `expr` must be of type `str`, sympy.Expr or callable, instead got {type(exprs[ii])}.")
            
    if not (all(isfunc) or all([not b for b in isfunc])):  # Must either be all expressions or all callables
        raise TypeError("`expr` cannot mix callable functions with expressions.")
        
    # Validate intervals is given in the correct format and the correct number given within the bounds.
    if is_listlike(intervals):
        if len(intervals) == 0:
            if len(exprs) == 1:  # if no intervals are given, set the interval to be the ending bound as a placeholder
                intervals = [bounds[-1]]
            else:
                raise ValueError(f"No intervals given for {len(exprs)} expressions.")
        elif len(intervals) != len(exprs) - 1:
            raise ValueError(f"There should be N-1 intervals for N expressions, instead received {len(intervals)} intervals for {len(exprs)} expressions.")
    else:
        intervals = [intervals]
        if len(intervals) != len(exprs) - 1:
            raise ValueError(f"There should be N-1 intervals for N expressions, instead received {len(intervals)} intervals for {len(exprs)} expressions.")
    
    # Ensure they are within the bounds given.
    if min(intervals) < bounds[0]:
        raise ValueError(f"The minimum interval value should be within the current bounds ({bounds[0]}, {bounds[1]}).")
    elif max(intervals) > bounds[1]:
        raise ValueError(f"The maximum interval value should be within the current bounds ({bounds[0]}, {bounds[1]}).")

    intervals = list(intervals)
    
    # Validate the substitution variable names
    symbolkeysArr = []
    if (var in list(subs)):
        raise ValueError(f"Independent variable '{var}' should not be in subs")
        
    if subs is None: subs = dict()
    if not any(isfunc):  # expr contains Sympy expressions
        # symbols = exprs[0].free_symbols  ## removing this line as probably not needed?
        for ii in range(len(exprs)):
            symbols = exprs[ii].free_symbols  # Get expression Symbols
            symbolkeys = [str(symb) for symb in symbols]  # convert these to strings, instead of sympy.Symbols
            
            if var is None:  # Assume independent variable is 'x', otherwise, assume the first symbol.
                if orientation == 'horizontal':
                    var = 'x'  # if 'x' in symbolkeys or len(symbolkeys)==0 else None
                else:  # first test for 'y' as an independent variable, then default to x.
                    if 'y' in symbolkeys:
                        var = 'y'
                    else:
                        var = 'x'  # if 'x' in symbolkeys or len(symbolkeys)==0 else symbolkeys[0]
                        
            # remove the independent variable from the symbols and append to list
            symbolkeysArr += symbolkeys
        
        # Check for missing substitution variables
        missing = set(symbolkeysArr) - set(subs) - set(var)
        if len(missing) > 0:
            raise TypeError(f"`expr` missing {len(missing)} required substitution variable{'s' if len(missing) > 1 else ''}: {list(missing)}")

        # Check for any substitution variables that are not in any expressions
        for key in list(subs):
            if not any([key in symb for symb in symbolkeysArr]):
                raise KeyError(f"Substitution variable '{key}' does not exist in any 'expr'")

    # The lengths of each substitute value list, len=1 if just a single value
    lens = [len(subs[key]) if (isinstance(subs[key], Iterable) and not isinstance(subs[key], str)) else 1 for key in list(subs)]
    if (permute is True):
        L = prod(lens)
        perms = array(meshgrid(*subs.values())).reshape(len(subs), -1)
        permsubs = {}
        for ii, key in enumerate(list(subs)):
            permsubs[key] = perms[ii]
        subsarr = simpler_dict_splicer(permsubs, L, [1] * L)
    else:
        L = max(lens) if len(lens) > 0 else 1
        subsarr = simpler_dict_splicer(subs, L, [1] * L)
    
    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    plot_par = {**plot_kw, **kwargs}

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    # Set default label behaiviour
    if label is None:
        label = False if isfunc else True

    # Create the labels
    print('label', label)
    if label is False:
        labellist = [None] * L
    elif label is True:  # Auto-generate labels
        if all(isfunc):
            labellist = []
            exprstr = "; ".join(str(expr.__name__) for expr in exprs)
            
            pars = []
            for expr in exprs:
                signature = inspect.signature(expr)  # get the signature of the function/lambda
                pars += [str(par) for par in list(signature.parameters.values())[1:] if str(par) not in pars]
            
            for ii in range(L):  # Make a label for each of sub values
                labellist.append(f"{exprstr} {{" + ", ".join(f"{par}={subsarr[ii][par]}" if par in list(subsarr[ii]) else f"{par}" for par in pars) + "}")
                
        elif subsarr == [{}]:
            labellist = ["; ".join(f"${latex(expr)}$" for expr in exprs) if uselatex else "; ".join(exprs)]
        else:
            labellist = []
            exprstr = "; ".join(f"${latex(expr)}$" for expr in exprs) if uselatex else "; ".join(exprs)
            for ii in range(L):  # Make a label for each of sub values
                if uselatex:
                    labellist.append(f"{exprstr} {{" + "; ".join([f"${key}$={subsarr[ii][key]}" for key in list(subsarr[ii])]) + "}")
                else:
                    labellist.append(f"{exprstr} {{" + "; ".join([f"{key}={subsarr[ii][key]}" for key in list(subsarr[ii])]) + "}")
    
    elif isinstance(label, str):  # A single string
        labellist = [label] * L
    
    else:
        if not is_listlike(label):
            raise TypeError(f"`label` of type {type(label)} is not recognised.")
            
        if (len(label) != L):
            raise TypeError(f"Number of labels ({len(label)}) does not match the number of curves ({L}).")
        
        labellist = label

    vararr = logspace(*log10(bounds), num=num) if xlog else linspace(*bounds, num=num)
    intervals.append(vararr[-1])
    curves = [None] * L

    for ii in range(L):
        lines = []
        start = 0
        for jj, expr in enumerate(exprs):
            end = array(abs(vararr - intervals[jj])).argmin()
            if isfunc[jj]:
                signature = inspect.signature(expr)
                pars = [str(par) for par in signature.parameters.values()]
                subs = {key: val for key, val in subsarr[ii].items() if key in pars}
                lines.append(list(zip(vararr[start: end + 1], expr(vararr[start: end + 1], **subs))))
            else:
                lamb = lambdify(var, expr.subs(subsarr[ii]), modules='numpy')
                if expr.subs(subsarr[ii]).is_constant():
                    func = lambda x: full_like(x, lamb(x))
                else:
                    func = lambda x: lamb(array(x))
                
                if orientation == 'horizontal':
                    lines.append(list(zip(vararr[start: end + 1], func(vararr[start: end + 1]))))
                else:
                    lines.append(list(zip(func(vararr[start: end + 1]), vararr[start: end + 1])))
                    
            start = end  # move to next interval
        
        curves[ii] = LineCollection(lines, label=labellist[ii], **plot_par[ii])
        ax.add_collection(curves[ii])

    _plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, gridpar, ax)
    
    if old_ax is not None:
        sca(old_ax)
    
    return(curves[0] if len(curves) == 1 else curves)


####################################
# 1D histogram and binned statistics
####################################
def hist(data, bin_type=None, bins=None, dens=True, cumul=None, scale=None, weights=None, hist_type=None, orientation='vertical',
         v=None, vstat=None, nmin=0, color=None, xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=None, title=None,
         xlabel=None, ylabel=None, label=None, lab_loc=0, ax=None, grid=None, plot_kw={}, output=None, **kwargs):
    """1D histogram function.

    The plotting is done with pyplot.plot(), so histograms are shown with interpolated curves instead
    of the more common stepwise curve. For this reason splotch.histstep is a better choice for
    small datasets.

    Parameters
    ----------
    data : array-like or list
        If list it is assumed that each elemement is array-like.
    bin_type : {'number','width','edges','equal'}, optional
        Defines how is understood the value given in bins: 'number' for the desired number of bins,
        'width' for the width of the bins, 'edges' for the edges of bins, and 'equal' for making
        bins with equal number of elements (or as close as possible). If not given it is inferred
        from the data type of bins: 'number' if int, 'width' if float and 'edges'if ndarray.
    bins : int, float, array-like or list, optional
        Gives the values for the bins, according to bin_type.
    dens :  bool or list, optional
        If false the histogram returns raw counts.
    cumul : bool or list, optional
        If true, produces a cumulative distribution instead of a histogram.
    scale : float or list, optional
        Scaling to be applied to the counts.
    weights : array-like or None, optional
        An array of weights with the same shape as data. For each value in data, it will only
        contribute its given weight towards the bin count (instead of 1). If dens is True, weights
        will also be normalised so that the integral over the density remains 1. Default: None
    hist_type : str, optional.
        Defines the type of histogram to be drawn. 'line' and 'step' produce lines, with the former
        drawing lines conecting the values of each bin positioned on their centre, and the latter
        drawing a stepwise line, with the edges of each step coinciding with the bin edges.
        'bar' produces a bar plot. All have filled version (i.e., 'linefilled'), which fills the
        space between the edges of the histogram and 0.
    orientation : {'vertical', 'horizontal'}, optional (default: 'vertical')
        If 'horizontal', histograms values will be plot along the x-axis instead of the
        default behaviour of plotting in the y-axis.
    v : array-like or list, optional
        If a valid argument is given in vstat, defines the value used for the binned statistics.
    vstat : str, function  or list, optional
        Must be or contain one of the valid str arguments for the statistics variable
        in scipy.stats.binned_statistic ('mean’, 'median’, 'count’, 'sum’, 'min’ or 'max’) or
        function(s) that takes a 1D array and outputs an integer or float.
    nmin : int, optional (default: 0)
        The minimum number of points required in a bin in order to be plotted.
    color : color or list-like of colors or None, default: None
        A color or sequence of strings to assign to each dataset given in `data`.
        If `hist_type` is 'bar' or any 'filled' variant, then `color` refers to the facecolor,
        whereas if `hist_type` is 'step', then color refers to the edgecolor. Note that `color`
        will be overwritten if `facecolor` or `edgecolor` are also specified. Only `color` is
        relevant when `hist_type` = 'line'. If `None` given, uses the next color in the
        current color cycler.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If true inverts the x-axis.
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    yinvert : bool or list, optional
        If true inverts the y-axis.
    xlog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
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
    grid : boolean, dict or None
        The grid behaviour, acts according to:
         * If boolean, the grid is set on (True) or off (False).
         * If dict, allows specific `.Line2D` properties of the grid to be set.
         * If None, use default grid parameters if this is the initial plotting call, otherwise do nothing.
    plot_par : dict, optional
        Passes the given dictionary as a kwarg to the plotting function.
    output : boolean, optional
        If True, returns the values and edges of the resulting histogram.

    Returns
    -------
    n : list
        List containing the arrays with the values for each histogram drawn. Only provided
        if output is True.
    bins_edges : list
        List containing the arrays with the bin edges for each of the histograms drawn.
        Only provided if output is True.
    """

    # Get the relevant axis
    if ax is None: ax = gca()
    old_ax = axes_handler(ax)

    gridpar = grid_handler(grid, ax)

    if not is_listlike(data):
        data = [data]
    elif isinstance(data, ndarray):
        if (len(shape(data)) == 1 and data.dtype is not dtype('O')):
            data = [data]

    L = len(data)

    if not is_listlike(bin_type):
        bin_type = [bin_type] * L
    if not is_listlike(bins) or (len(shape(bin_type)) == 1):
        if bins is not None:
            bins = [bins] * L
        else:
            bins = [int((len(d))**0.4) for d in data]
    if not is_listlike(weights) or (len(shape(weights)) == 1):
        weights = [weights] * L
    if not is_listlike(dens):
        dens = [dens] * L
    if not is_listlike(cumul):
        cumul = [cumul] * L
    if not is_listlike(scale):
        scale = [scale] * L
    if not is_listlike(v) or (len(shape(v)) == 1):
        v = [v] * L
    if not is_listlike(nmin) or (len(shape(nmin)) == 1):
        nmin = [nmin] * L
    if not is_listlike(vstat):
        vstat = [vstat] * L
    if not is_listlike(label):
        label = [label for i in range(L)]
    if is_number(xlim):
        xlim = [nanmean(data) - xlim * nanstd(data), nanmean(data) + xlim * nanstd(data)]
    if ylog is None:
        ylog = Params.hist1D_yaxis_log
    if hist_type is None:
        hist_type = Params.hist1D_histtype
    if output is None:
        output = Params.hist1D_output

    if not is_listlike(hist_type):
        hist_type = [hist_type] * L

    xlim = lims_handler(xlim, ax)
    ylim = lims_handler(ylim, ax)

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    plot_par = {**plot_kw, **kwargs}

    # Check if width is given as a kwarg
    if 'width' in plot_par.keys():
        warn('Received kwarg width, this will be ignored in the histogram', UserWarning)
        if hist_type != 'bar':
            del plot_par['width']

    if color is None:  # check color as per matplotlib.pyplot.hist()
        if 'c' in plot_par.keys():
            color = plot_par.pop('c')
        else:
            color = [ax._get_lines.get_next_color() for i in range(L)]

    color = to_rgba_array(color)
    if len(color) not in (1, L):  # either length of 1 or length of number of datasets.
        raise ValueError(f"The `color`/`c` keyword argument should contain one color or one color per dataset, "
                         f"but {L} datasets and {len(color)} colors given.")

    if len(color) == 1:  # duplicate color for all datasets if only one given
        color = [color] * L

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    plot_type = {'line': plt_plot, 'linefilled': fill_between if orientation == 'vertical' else fill_betweenx,
                 'step': stairs, 'stepfilled': stairs,
                 'bar': bar if orientation == 'vertical' else barh, 'barfilled': bar if orientation == 'vertical' else barh}
    hist_centre = {'line': True, 'linefilled': True,
                   'step': False, 'stepfilled': False,
                   'bar': False, 'barfilled': False}

    bin_edges = []
    n_return = []

    for i in range(L):
        temp_data, bins_hist, bins_plot = bin_axis(data[i], bin_type[i], bins[i], log=xlog, plot_centre=hist_centre[hist_type[i]])
        n_check = histogram(temp_data, bins=bins_hist, density=False)[0]
        n_check = n_check >= nmin[i]

        if vstat[i]:
            temp_y = binned_statistic(temp_data, v[i], statistic=vstat[i], bins=bins_hist)[0]
        else:
            temp_y = histogram(temp_data, bins=bins_hist, density=False if scale[i] else dens[i], weights=weights[i])[0]

        if cumul[i]:
            temp_y = nancumsum(temp_y)
            if dens[i]:
                temp_y = temp_y.astype('float') / nanmax(temp_y)

        if scale[i]:
            temp_y = temp_y.astype('float') / scale[i]
            if dens[i]:
                temp_y /= bins_hist[1:] - bins_hist[:-1]
        if ylog:
            temp_y = where(temp_y == 0, nan, temp_y)

        temp_y = where(n_check, temp_y, nan)
        y = temp_y

        if hist_type[i] in ['step', 'line']:
            for par in ['edgecolor', 'ec']:  # edgecolor supercedes color for step and line funcs.
                if (par in plot_par[i].keys()):
                    color[i] = to_rgba_array(plot_par[i].pop(par))
            for par in ['facecolor', 'fc']:  # facecolor is superfluous for step and line funcs.
                if (par in plot_par[i].keys()):
                    del plot_par[i][par]

        if hist_type[i].startswith('step'):
            plot_par[i]['fill'] = True if hist_type[i] == 'stepfilled' else False
            plot_par[i]['orientation'] = orientation

        if hist_type[i].startswith('bar'):
            plot_par[i]['fill'] = True if hist_type[i] == 'barfilled' else False

            if orientation == 'vertical':
                plot_par[i]['width'] = diff(bins_plot)
            else:
                plot_par[i]['height'] = diff(bins_plot)

            bins_plot = (bins_plot[0:-1] + bins_plot[1:]) / 2

            if ('edgecolor' not in plot_par[i].keys() and 'ec' not in plot_par[i].keys()):  # could also be 'ec'
                if (color[i] is not None):
                    plot_par[i]['edgecolor'] = color[i]
                else:
                    cycler = ax._get_lines.prop_cycler
                    plot_par[i]['edgecolor'] = next(cycler)['color']
        
        if (orientation == 'vertical' and hist_type[i] not in ('step', 'stepfilled')) or hist_type[i] in ('bar', 'barfilled', 'linefilled'):
            plot_type[hist_type[i]](bins_plot, y, color=color[i], label=label[i], **plot_par[i])

        else:
            plot_type[hist_type[i]](y, bins_plot, color=color[i], label=label[i], **plot_par[i])

        bin_edges.append(bins_plot)
        n_return.append(temp_y)

    if any(label) and Params.legend_auto is True:
        legend(loc=lab_loc)

    if ylim == 'auto' and orientation == 'vertical':  # Adjust ylims if None given.
        if not ylog and all([val is None for val in v]):  # These automatic limits do not apply when ylog=True or statistics are used.
            ylim = [0, max(nanmax(y) * (1 + rcParams['axes.ymargin']), gca().get_ylim()[1])]

    if xlim == 'auto' and orientation == 'horizontal':  # Adjust xlims if 'auto'.
        if not xlog and all([val is None for val in v]):  # These automatic limits do not apply when ylog=True or statistics are used.
            xlim = [0, max(nanmax(y) * (1 + rcParams['axes.xmargin']), gca().get_xlim()[1])]

    _plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, gridpar, ax)

    if old_ax is not ax:
        sca(old_ax)

    if len(n_return) == 1:
        n_return = n_return[0]

    if len(bin_edges) == 1:
        bin_edges = bin_edges[0]

    if output:
        return(n_return, bin_edges)


####################################
# Standard plots
####################################
def plot(x, y=None, xlim=None, ylim=None, xinvert=False, yinvert=False, xlog=False, ylog=False, title=None, xlabel=None,
         ylabel=None, label=None, lab_loc=0, ax=None, grid=None, plot_kw={}, **kwargs):
    """Base plotting function.

    This is a wrapper for pyplot.plot().

    Parameters
    ----------
    x : array-like or list
        If list it is assumed that each elemement is array-like. If y is not given, the given
        values pass to y and a numpy array is generated with numpy.arange() for the x values.
    y : array-like or list, optional
        If list it is assumed that each elemement is array-like.
    xlim : tuple-like, optional
        Defines the limits of the x-axis, it must contain two elements (lower and higer limits).
    ylim : tuple-like, optional
        Defines the limits of the y-axis, it must contain two elements (lower and higer limits).
    xinvert : bool or list, optional
        If true inverts the x-axis.
    yinvert : bool or list, optional
        If true inverts the y-axis.
    xlog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    ylog : bool or list, optional
        If True the scale of the x-axis is logarithmic.
    title : str, optional
        Sets the title of the plot
    xlabel : str, optional
        Sets the label of the x-axis.
    ylabel : str, optional
        Sets the label of the y-axis.
    label : str, optional
        Sets the legend for the plot.
    lab_loc : int, optional
        Defines the position of the legend
    ax : pyplot.Axes, optional
        Use the given axes to make the plot, defaults to the current axes.
    grid : boolean, dict or None
        The grid behaviour, acts according to:
         * If boolean, the grid is set on (True) or off (False).
         * If dict, allows specific `.Line2D` properties of the grid to be set.
         * If None, use default grid parameters if this is the initial plotting call, otherwise do nothing.
    plot_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are
        Line2D properties.
    **kwargs: Line2D properties, optional
        kwargs are used to specify matplotlib specific properties such as linecolor, linewidth,
        antialiasing, etc. A list of available `Line2D` properties can be found here:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    Returns
    -------
    lines
        A list of Line2D objects representing the plotted data.
    """

    # Get the relevant axis
    if ax is None: ax = gca()
    old_ax = axes_handler(ax)

    gridpar = grid_handler(grid, ax)
    xlim = lims_handler(xlim, ax)
    ylim = lims_handler(ylim, ax)

    if not is_listlike(x) or len(shape(x)) == 1:
        x = [x]

    L = len(x)
    
    if y is None:
        y = x
        x = [arange(len(x[i])) for i in range(L)]
    else:
        if not is_listlike(y) or len(shape(y)) == 1:
            y = [y]

    if not is_listlike(label):
        label = [label for i in range(L)]

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    plot_par = {**plot_kw, **kwargs}

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    plot_par = dict_splicer(plot_par, L, [1] * L)

    lines = []  # Initialising list which contains each line
    for i in range(L):
        lines += ax.plot(x[i], y[i], label=label[i], **plot_par[i])

    if any(label):
        legend(loc=lab_loc)

    _plot_finalizer(xlog, ylog, xlim, ylim, title, xlabel, ylabel, xinvert, yinvert, gridpar, ax)

    if old_ax is not ax:
        sca(old_ax)

    return (lines[0] if len(lines) == 1 else lines)
