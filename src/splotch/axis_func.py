def adjust_text(which=['x', 'y'], ax=None, text_kw={}, **kwargs):
    """ Adjusts text instances.

    Function which allows the user to adjust texts of one or many of either x/y-axis labels,
    title, legend, colorbars, etc.

    Parameters
    ----------
    which : str, array-like or list
        Which ``Text`` instance(s) to apply the text properties to. Valid arguments can be one or many of the following

        ==============================  =================================
         Arguments                      Description
        ==============================  =================================
        'x','xlabel'                     x-axis label
        'y', 'ylabel'                    y-axis label
        'k', 'tick'                      Tick labels
        't', 'title'                     Title
        's', 'suptitle'                  Sup. title
        'l', 'legend'                    Legend text
        'c', 'colorbar'                  Color bar
        'T', 'text'                      Text objects
        'a', 'all'                       All instances of all the above
        ==============================  =================================

    ax : pyplot.Axes or list, optional
        Use the given axes to adjust text, defaults to the current axis. If a list of Axis instances given, the text properties is applied to each.

    text_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib ``Text`` instance. It is recommended that text keyword arguments be given as **kwargs.

    **kwargs : Text instance properties
        kwargs are used to specify properties of Text instances. A list of valid Text kwargs can be found in the matplotlib
        `Text <https://matplotlib.org/stable/api/text_api.html>`_ documentation.

    """

    from matplotlib.text import Text
    from matplotlib.pyplot import gca
    from numpy import array, shape, max as np_max, argmax, append
    from .base_func import axes_handler, dict_splicer, plot_finalizer

    try:
        _ = (it for it in ax)
    except TypeError:
        if (ax is None):
            ax = gca()
        ax = [ax]

    # Validate `which` value(s)
    whichRef = ['x', 'y', 't', 's', 'k', 'l', 'c', 'T', 'a',
                'xlabel', 'ylabel', 'title', 'suptitle', 'ticks', 'legend', 'colorbar', 'text', 'all']

    try:  # check if iterable
        _ = (i for i in which)
        if (isinstance(which, str)):  # If string instance
            if which in whichRef[9:]:  # Set to a single-item list if one of the long references
                which = [which]
            else:  # Else, check this is not a combination of shortened references i.e., 'xyk'
                # Check condition: All of the characters are recognised letters and there are no duplicates
                if (all(w in whichRef[:9] for w in list(which)) and len(which) <= len(set(which))):
                    which = list(which)  # set to a list of each of the letters
                else:
                    which = [which]

    except (TypeError):
        which = [which]

    for w in which:
        try:
            whichIndex = whichRef.index(w)
            whichComp = (whichIndex + len(whichRef) // 2) % len(whichRef)  # get the complimenting short/long version
            if (whichRef[whichComp] in which):
                raise TypeError("adjust_text() received equivalent values for 'which': '{0}' and '{1}'.".format(whichRef[whichIndex], whichRef[whichComp]))
        except (ValueError):
            if (not isinstance(w, Text)):
                raise TypeError("adjust_text() received invalid value for 'which' ('{0}'). Must be one of: {1}".format(w, ', '.join(whichRef)))

    L = len(which)

    # Combine the `explicit` plot_kw dictionary with the `implicit` **kwargs dictionary
    # plot_par={**plot_kw, **kwargs} # For Python > 3.5
    textpar = text_kw.copy()
    textpar.update(kwargs)

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    textpar = dict_splicer(textpar, L, [1] * L)

    for a in array(ax).flatten():
        if (a is None or a.axes is None):
            continue  # ignore empty subplots

        for ii, lab in enumerate(which):
            if (lab in ['x', 'xlabel']):
                texts = [a.xaxis.label]
            elif (lab in ['y', 'ylabel']):
                texts = [a.yaxis.label]
            elif (lab in ['t', 'title']):
                texts = [a.title]
            elif (lab in ['s', 'suptitle']):  # Not implemented
                texts = [a.title]
            elif (lab in ['k', 'ticks']):
                texts = append(a.get_yticklabels(), a.get_xticklabels())
            elif (lab in ['l', 'legend']):
                lgnd = a.get_legend()
                if (lgnd is not None):
                    texts = lgnd.get_texts()  # Get list of text objects in legend (if it exists)
            elif (lab in ['c', 'colorbar']):
                # Get the axis with the largest ratio between width or height
                caxInd = argmax([np_max([c.get_position().width / c.get_position().height, c.get_position().height / c.get_position().width]) for c in a.figure.get_axes()])
                if (a.figure.get_axes()[caxInd].get_position().height > a.figure.get_axes()[caxInd].get_position().width):
                    texts = [a.figure.get_axes()[caxInd].yaxis.label]
                else:
                    texts = [a.figure.get_axes()[caxInd].xaxis.label]
            elif (lab in ['T', 'text']):
                texts = [child for child in a.get_children()[:-4] if isinstance(child, Text)]  # -4 to avoid grabbing Title, Subtitle. etc. Text instances
            elif (isinstance(lab, Text)):
                texts = [lab]
            elif (lab in ['a', 'all']):
                texts = [a.xaxis.label, a.yaxis.label, a.title, *a.get_xticklabels(), *a.get_yticklabels()]

                lgnd = a.get_legend()
                if (lgnd is not None): texts = texts + lgnd.get_texts()

                caxInd = argmax([np_max([c.get_position().width / c.get_position().height, c.get_position().height / c.get_position().width]) for c in a.figure.get_axes()])
                if (a.figure.get_axes()[caxInd].get_position().height > a.figure.get_axes()[caxInd].get_position().width):
                    texts.append(a.figure.get_axes()[caxInd].yaxis.label)
                else:
                    texts.append(a.figure.get_axes()[caxInd].xaxis.label)

                texts = texts + [child for child in a.get_children()[:-4] if isinstance(child, Text)]

            for t in texts:  # Actually apply the font changes
                t.set(**textpar[ii])

    return(None)


def colorbar(mappable=None, ax=None, label='', orientation='vertical', loc=1, transform=None,
             inset=False, aspect=0.05, width=None, height=None, pad=0.05, ticks=None, bar_kw={}, **kwargs):
    """Colorbar function

    This function will produce a colorbar for any currently plotted mappable on a given axis.

    Parameters
    ----------
    mappable : ScalarMappable, optional
        A list or individual matplotlib.cm.ScalarMappable (i.e., Image, ContourSet, etc.) to be described
        by this colorbar(s). If no mappable given, any currently present mappables will be searched for
        in each axis given.
    ax : pyplot.Axes, optional
        Use the axis specified to draw the colorbars onto. If multiple axes given, the number of
        mappable objects provides must be one or equal to the number of axes.
        Defaults to the current axis.
    label : str, optional
        The label to be given to the colorbar.
    orientation : optional
        The orientation of the colorbar, permitted values are: 'vertical' or 'horizontal'.
        Orientation is necessary to decide which spine of the colorbar axis to place labels.
        Default: 'vertical'.
    loc : int or tuple-like, optional
        Specifies the location of the colorbar. Can either be the string or integer case from the list below.

            ===============   =============
            Location String   Location Code
            ===============   =============
            \'upper right\'     1
            \'upper left\'      2
            \'lower left\'      3
            \'lower right\'     4
            \'center right\'    5
            \'upper center\'    6
            \'center left\'     7
            \'lower center\'    8
            \'center\'          9
            ===============   =============

    transform : matplotlib.transforms.Transform instance, optional
        The transformation instance to be used for colorbar location if loc is tuple-like. For example,
        using ax.transAxes() will specify the colorbar location in the coordinates of the axis; (0,0)
        is bottom-left and (1,1) is top-right. Default: ax.transAxes.
    inset : boolean, optional
        Whether to inset the colorbar within the inside of the axis. If loc not tuple-like, this
        will add padding to both sides of the colorbar to avoid colliding with the axis spine.
        Default: False.
    aspect : float, optional
        The aspect ratio of the colorbar always taken as the ratio of the long-side to the
        short-side. Default: 0.05
    width, height : float, optional
        The width and height of the colorbar in the coordinate system of specified transform,
        which by default is in the coordinates of the axis.
    pad : float, optional
        The padding given to the colorbar axis offset from the margin of the axis. If loc is
        specified to an edge which has an axis label/ticks, padding will be added from the edge
        of the label.
    ticks : None, list-like or Locator() object, optional
        Specifies the locations of ticks on th colorbar axis. If None, ticks are determined
        automatically from the input.
    bar_kw : dict, optional
        Passes the given dictionary as a kwarg to the plotting function. Valid kwargs are colorbar
        properties.
    **kwargs : Colorbar properties, optional
        Keyword arguments are used to specify matplotlib.pyplot.colorbar specific properties such as
        extend, spacing, format, drawedges, etc. The list of available properties can be found in the
        matplotlib `Colorbar <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html>`_ documentation.

    Returns
    -------
    cbar
        The colorbar object.
    """

    from splotch.base_func import axes_handler, dict_splicer
    from matplotlib.pyplot import gca, colorbar
    from matplotlib.contour import ContourSet
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib.ticker as tckr
    from numpy import ndarray

    # Validate axis input
    if ax is not None:
        try:  # check if iterable
            _ = (i for i in ax)
            # old_axes = axes_handler(ax[0])
            axes = ax
        except (TypeError):
            # old_axes = axes_handler(ax)
            axes = [ax]
    else:
        axes = [gca()]
        # old_axes = axes[0]

    locRef = {'upper right': 1, 'top right': 1,
              'upper left': 2, 'top left': 2,
              'lower left': 3, 'bottom left': 3,
              'lower right': 4, 'bottom right': 4,
              'right': 5, 'center right': 5,
              'upper center': 6, 'top center': 6,
              'left': 7, 'center left': 7,
              'lower center': 8, 'bottom center': 8,
              'center': 9}

    if isinstance(loc, str):
        loc = loc.replace('centre', 'center')  # allows both spellings of centre/center
        if loc in locRef.keys():
            loc = locRef[loc]
        else:
            raise ValueError(f"loc {loc} not recognised.")
    elif isinstance(loc, int):
        if loc < 0 or loc > 9:
            raise ValueError(f"loc {loc} must be between 0 - 9.")
    else:
        raise NotImplementedError("loc must be specified as integer or string. Providing loc as a tuple-like as colorbar anchor position is not yet implemented.")

    ### Define the positions of preset colorbars
    labpad = 0.125  # The padding added for labels
    ins = 1 if inset is True else 0  # Convert inset boolean to a binary multiplier

    # Validate aspect value
    if (aspect <= 0 or aspect > 1):
        raise ValueError("Value for aspect must be strictly positive and less than or equal to 1 (i.e. 0 < aspect <= 1)")

    # Get width/height values of colorbar
    if (width is None and height is None):
        width, height = (aspect, 1.0 - ins * 2 * pad) if orientation == 'vertical' else (1.0 - ins * 2 * pad, aspect)
    else:
        if height is None:
            height = width / aspect if orientation == 'vertical' else aspect * width
        elif width is None:
            width = aspect * height if orientation == 'vertical' else height / aspect

    # Vertically-oriented colorbars
    vertPositions = {1: (1 + pad - ins * (2 * pad + width), 1 - height - ins * pad, width, height),  # upper right
                     2: (0 - width - labpad - pad + ins * (labpad + 2 * pad + width), 1 - height - ins * pad, width, height),  # upper left
                     3: (0 - width - labpad - pad + ins * (labpad + 2 * pad + width), 0 + ins * pad, width, height),  # lower left
                     4: (1 + pad - ins * (2 * pad + width), 0 + ins * pad, width, height),  # lower right
                     5: (1 + pad - ins * (2 * pad + width), 0.5 * (1 - height), width, height),  # center right
                     6: (0.5 * (1 - width), 1 + pad - ins * (2 * pad + height), width, height),  # upper center
                     7: (0 - width - labpad - pad + ins * (labpad + 2 * pad + width), 0.5 * (1 - height), width, height),  # center left
                     8: (0.5 * (1 - width), 0 - height - labpad - pad + ins * (labpad + 2 * pad + height), width, height),  # lower center
                     9: (0.5 * (1 - width), 0.5 * (1 - height), width, height)}  # center

    # horizontally - oriented colorbars
    horPositions = {1: (1 - width - ins * pad, 1 + pad - ins * (2 * pad + height), width, height),
                    2: (0 + ins * pad, 1 + pad - ins * (2 * pad + height), width, height),  # upper left
                    3: (0 + ins * pad, 0 - height - labpad - pad + ins * (height + labpad + 2 * pad), width, height),  # lower left
                    4: (1 - width - ins * pad, 0 - height - labpad - pad + ins * (height + labpad + 2 * pad), width, height),  # lower right
                    5: (1 + pad - ins * (2 * pad + width), 0.5 * (1 - height), width, height),  # center right
                    6: (0.5 * (1 - width), 1 + pad - ins * (2 * pad + height), width, height),  # upper center
                    7: (0 - width - labpad - pad + ins * (width + labpad + 2 * pad), 0.5 * (1 - height), width, height),  # center left
                    8: (0.5 * (1 - width), 0 - height - labpad - pad + ins * (labpad + 2 * pad + height), width, height),  # lower center
                    9: (0.5 * (1 - width), 0.5 * (1 - height), width, height)}  # center

    # Combine the `explicit` bar_kw dictionary with the `implicit` **kwargs dictionary
    # bar_par={**bar_kw, **kwargs} # For Python > 3.5
    bar_par = bar_kw.copy()
    bar_par.update(kwargs)

    # Create 'L' number of plot kwarg dictionaries to parse into each plot call
    bar_par = dict_splicer(bar_par, len(axes), [1] * len(axes))

    cbars = []  # Initiate empty list of colorbars for output
    for ii, ax in enumerate(axes):
        if mappable is None:
            # Try to automatically find mapable as children of the axis
            mappables = [child for child in ax.get_children() if hasattr(child, 'autoscale_None')]
            if mappables != []:  # A potential mappable was detected
                mapper = mappables[0] if isinstance(mappables, (tuple, list, ndarray)) else mappables
            else:
                if ax.figure._gci() is not None:
                    mapper = ax.figure._gci()
                else:
                    raise ValueError("No potential mappables found in figure.")
        else:
            # check if list-like
            if isinstance(mappable, (list, tuple, ndarray)):
                if (len(mappable) == len(axes)):
                    mapper = mappable[ii]
                elif (len(mappable) == 1):
                    mapper = mappable[0]
                else:
                    raise ValueError("Number of mappables given must be either 1 or equal to the number of axes specified.")
            else:
                mapper = mappable

            if isinstance(mapper, ContourSet) and not mapper.filled:
                raise NotImplementedError("Cannot create colorbars for unfilled contours.")

        cax = inset_axes(ax, width='100%', height='100%',
                         bbox_to_anchor=vertPositions[loc] if orientation == 'vertical' else horPositions[loc],
                         bbox_transform=transform if transform is not None else ax.transAxes,
                         borderpad=0)

        cbar = colorbar(mapper, cax=cax, orientation=orientation, ticks=ticks, **bar_par[ii])
        cbar.ax.yaxis.set_major_formatter = tckr.ScalarFormatter()

        # Orient tick axes correctly
        if (orientation == 'horizontal'):
            cbar.ax.set_xlabel(label)

            flip = True if loc in [3, 4, 8] else False  # Flip the labels if colorbar on bottom edge
            if (inset is True) and loc not in [5, 7, 9]:
                flip = not flip  # Reverse flipping in the case of a inset colorbar

            if (flip is True):
                cbar.ax.xaxis.set_label_position('bottom')
                cbar.ax.xaxis.tick_bottom()
            else:
                cbar.ax.xaxis.set_label_position('top')
                cbar.ax.xaxis.tick_top()

            cbar.ax.yaxis.set_ticks([], minor=True)

        else:
            cbar.ax.set_ylabel(label)

            flip = True if loc in [2, 3, 7] else False  # Flip the labels if colorbar on left edge
            if (inset is True) and loc not in [6, 8, 9]:
                flip = not flip  # Reverse flipping in the case of a inset colorbar

            if (flip is True):
                cbar.ax.yaxis.set_label_position('left')
                cbar.ax.yaxis.tick_left()

            cbar.ax.xaxis.set_ticks([], minor=True)

        cbars.append(cbar)

    return (cbars[0] if len(cbars) == 1 else cbars)


def cornerplot(data,columns=None,pair_type='contour',nsamples=None,sample_type='rand',labels=None,histlabel=None,
                fig=None,figsize=None,wspace=0.0,hspace=0.0,squeeze=False,
                hist_kw={},contour_kw={},scatter_kw={},hist2D_kw={},axes_kw={},_debug_=False,**kwargs):
    """ Creates a corner plot figure and subplots
    
    This function accepts columns of data representing multiple parameters in which each combination
    will be paired will each other to form the off-diagonals of a corner plot. The diagonals show a
    1D histogram representation of each individual parameter.
    
    Parameters
    ----------
    data : array-like
        The input data frame with each parameter represented by an individual column, the zeroth axis should be the list of samples and the next axis should be the number of dimensions. Accepted data types are: pandas.DataFrame, pandas.Series, numpy.ndarray, astropy.table.Table.

    columns : array-like
        The column labels (or indices) that specify which columns within 'data' to use. if none specified, every column in 'data' with a numeric datatype will be used. To group columns together, columns can be given as a list of lists, e.g.:
            columns=[['A1', 'A2', ...], ['B1', 'B2', ...]]

        which will create pairs for all parameters in each sublist. Having multiple groups will result in multiple histograms and paired plots per axis. As such, pair_type='hist2D' is not compatible when columns includes multiple groups.
    pair_type : str, optional
        The plotting type for the off-diagonal plots,
        can be one of:'contour' | 'scatter' | 'hist2D' which correspond to contour plots, scatter
        plots and 2D histograms. A dictionary of parameters can be parsed to 'pair_kw' to control
        the styling of each.
    nsamples : float, optional
        Specifies the number of samples kept out of the total number of samples. This is useful to
        speed up plotting if not all of the samples need to be shown. Default: 1 (i.e. no subsampling).
        Default: the full number samples in ``data``.
    sample_type : str, optional, Default='end'
        Sets the method used to thin the full set of samples, can be one of the following:
            * 'end' : Takes the last ``nsamples`` samples (Useful for MCMC posteriors where the end is often better).
            * 'rand' : Randomly selects a set of samples. (Only use if confident all posterior chains are stationary).
            * 'thin' : Evenly select every m samples until a total of ``nsamples`` are kept.

    labels : array-like (str), optional
        A list of axis labels to be assigned to each parameter. Must be of the same length or longer
        than the number of columns (or column groups).
    histlabel : str, optional
        The y-axis label for each of the 1D histograms on the diagonal, if None given, then no
        labels or ticks will be drawn. Labels and ticks will be displayed on the opposite side to
        the labels for the pair plots, Default: None.
    figsize : 2-tuple of floats, default: rcParams["figure.figsize"] * (len(data), len(data))
        The dimensions of the figure (width, height) in inches. If not specified, the default is to
        scale the default rcParams figure.figsize by the number of rows or columns.
    wspace / hspace : float, optional
        The horzontal/vertical spacing between figure subplots, expressed as a fraction of the
        subplot width/height.
    squeeze : bool, optional, default: True
        As per matplotlib's usage, the following applies:
            * If True, extra dimensions are squeezed out from the returned array of Axes: For NxN, subplots with N>1 are returned as a 2D array or as a scalar otherwise. 2D arrays fill out subplots from top-to-bottom and left-to-right.
            * If False (Default), no squeezing at all is done: the returned Axes object is always a 2D array containing Axes instances, even if it ends up being 1x1. For cornerplots with only one set of diagonals, empty axes will be filled by None.

    contour_kw, scatter_kw, hist2D_kw : dict, optional
        Dictionary of keyword arguments to be parsed into the pair plotting functions. The particular
        keyword dictionary used is dependent on the 'pair_type' chosen. If 'columns' contains
        multiple groups, kw arguments in lists are assumed to correspond to each group.
    hist_kw : dict, optional
        Dictionary of keyword arguments to be parsed into the 1D histograms plots on the diagonals.
    \*\*kwargs : Subplot instance properties
        kwargs are used to specify properties of ``subplots`` instances
        A list of valid ``axis`` kwargs can be found here:
        [https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes](https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes "Matplotlib.axes.Axes")    
    """
    from splotch.base_func import is_numeric, dict_splicer
    from splotch.plots_2d import contourp, scatter, hist2D
    from numpy import shape, reshape, full, ndarray, array, arange
    from numpy.random import choice
    from pandas import DataFrame, Series, RangeIndex
    
    from matplotlib.pyplot import figure
    from matplotlib.figure import Figure
    from matplotlib import rcParams
    from matplotlib.gridspec import GridSpec
    
    import warnings
    
    # Try to import astropy
    try:
        from astropy.table import Table
        hasAstropy=True
    except ImportError:
        hasAstropy=False
    
    nGroups=shape(columns)[0] if len(shape(columns)) > 1 else 1    
    if (_debug_ == True): print(f"Groups: {nGroups}")
            
    # Get number of parameters and dimensions
    dims=shape(data)
    
    # Validate plot parameters
    if (isinstance(pair_type,(list, ndarray, tuple))):
        if len(pair_type) > 2:
            raise ValueError("Too many pair types given.")
    elif (isinstance(pair_type, str)):
        pair_type=[pair_type]
    else:
        raise ValueError(f"Object type '{type(pair_type)}' not recognised for 'pair_type' argument. Must be either string of list-like.")
    for pair in pair_type:
        if (pair not in ['contour','scatter','hist2D']):
            raise ValueError(f"Pair type '{pair}' not a valid argument. Must be one of {{'contour'|'scatter'|'hist2D'}}.")
        if (pair == 'hist2D' and nGroups > 1):
            raise ValueError(f"Cannot overplot groups of 2D histograms. Choose alternative 'pair_type'.")
    
    # Sample the data
    if (nsamples == None):
        nsamples=dims[0] # By default, use all samples.
    else:
        if (nsamples > dims[0]):
            raise ValueError(f"Number of samples ({nsamples}) is greater than the total number of samples in data ({dims[0]})")
        elif (nsamples <= 0):
            raise ValueError(f"Number of samples ({nsamples}) must be positive and non-zero.")
    if (sample_type.lower() == 'end'):
        samps=arange(dims[0]-nsamples,dims[0])
    elif (sample_type.lower() == 'rand'):
        samps=choice(arange(0,dims[0]),size=nsamples,replace=False)
    elif (sample_type.lower() == 'thin'):
        step=(dims[0]-1)//(nsamples-1) if nsamples > 1 else dims[0] # Get the largest possible step size given the sample size
        offset=(dims[0] - step*(nsamples-1) - 1) // 2 # Offset of the initial point so that the points are centered
        samps=array([offset + kk*step for kk in range(nsamples)])
    else:
        raise ValueError(f"Sample type '{sample_type}' not recognised.")
    
    if (_debug_ == True): print(f"Number of samples: {len(samps)}")
    
    # Validate input data
    if isinstance(data, DataFrame) or (hasAstropy and isinstance(data, Table)):
        if (hasAstropy and isinstance(data, Table)): # Convert astropy.Table to pandas DataFrame if required
            data=data.to_pandas()
        if (columns == None): # no specific columns given
            cols=[c for c in data.columns if is_numeric(data[c])] # only include columns that are strictly numeric
            columns=cols
        else:
            cols=[]; nonNumeric=[]
            for kk, c in enumerate(array(columns).flatten()):
                if (is_numeric(data[c])):
                    cols.append(c)
                else:
                    nonNumeric.append(c)
                    del labels[kk]
            
            if (len(nonNumeric) > 0): # Warn if any columns were not numeric
                warnings.warn("Data type of column(s) '{0}' not numeric, ignoring column(s).".format(','.join(nonNumeric))) 
            
        cols=reshape(cols, shape(columns)).T
        
        # check that at least one numeric column was found
        if (len(cols) == 0): raise ValueError("No numeric columns found in data.")
        
    elif isinstance(data, Series):
        data=DataFrame(data)
        cols=list(data.columns)
    elif isinstance(data, ndarray):
        data=DataFrame(data) # recast numpy array into pandas DataFrame
        
        if (columns == None): # no specific columns given
            cols=[c for c in data.columns if is_numeric(data[c])] # only include columns that are strictly numeric
        else:
            cols=[]; nonNumeric=[]
            for kk, c in enumerate(array(columns).flatten()):
                if (is_numeric(data[c])):
                    cols.append(c)
                else:
                    nonNumeric.append(c)
                    del labels[kk]
            
            if (len(nonNumeric) > 0): # Warn if any columns were not numeric
                warnings.warn("Data type of column(s) '{0}' not numeric, ignoring column(s).".format(','.join(nonNumeric))) 
        
        cols=reshape(array(data[cols].columns), shape(cols)).T
    else:
        try: # Attempt to cast input data of unknown type into pandas DataFrame
            data=DataFrame(data)
        except (ValueError):
            raise ValueError(f"'data' must be pandas DataFrame/Series, np.ndarray or astropy.Table object, not: {type(data)}")
    
    
    # Get the number of parameters to create axes for
    npar=cols.size//nGroups if len(dims) > 1 else 1 # second axis defines the dimensions of parameters
    if (_debug_ == True): print(f"\nDimensions: {dims}")
    if (_debug_ == True): print(f"\nAxes: {npar}")
    
    # Assign labels if none given
    if (labels == None): # auto-generate labels from columns if available
        if not isinstance(data.columns, RangeIndex):
            labels=cols if nGroups<=1 else [c[0] for c in cols]
    else:
        if (len(labels) < npar):
            raise ValueError(f"Not enough labels given to match dimensions of data ({dims})")
    
    if (_debug_ == True): print(f"\nColumns:\n {cols}")
    
    #if (figsize==None): # Auto scale default figure size to num cols/rows.
    #    figsize=(rcParams["figure.figsize"][0]*npar, rcParams["figure.figsize"][0]*npar)
    
    axes=full(shape=(npar,npar),fill_value=None)
    if (fig is None): # Initialise figure and empty array
        fig=figure(figsize=figsize,**kwargs)
    else: # figure object was given
        if (fig.axes != []): # figure does contain axis objects
            try:
                if (len(fig.axes) != sum([kk for kk in range(1,npar+1)])):
                    raise ValueError(f"Number of axes in `fig` ({len(fig.axes)}) cannot be mapped to the number of axes required ({sum([kk for kk in range(1,npar+1)])}) for the given parameters.")
                else: # Reshape the new axes matrix
                    cnt=0
                    for ii in range(npar-1, -1, -1): # reshape axes list into appropriate matrix
                        for jj in range(ii if len(pair_type) == 1 else npar-1, -1, -1):
                            axes[ii,jj]=fig.axes[cnt]
                            cnt += 1
            except (TypeError): # axes in fig was not a list-like object
                raise TypeError("Figure axes has incorrect type, should be a list of 'matplotlib.axes.Axes' objects.")
    
    if (all([ax is None for ax in axes.flatten()])):
        # Set up the GridSpec object
        gs=GridSpec(ncols=npar,nrows=npar,wspace=wspace,hspace=hspace)
        for ii in range(npar-1, -1, -1):
            for jj in range(ii if len(pair_type) == 1 else npar-1, -1, -1):
                axes_kw['sharex']=None if ii==npar-1 else axes[npar-1,jj]
                if ii==jj: # Do not share y-axes of 1D histograms
                    axes_kw["sharey"]=None
                elif ii==npar-1: # 
                    axes_kw["sharey"]=None if jj==npar-2 else axes[ii,npar-2]
                else:
                    axes_kw["sharey"]=None if jj==npar-1 else axes[ii,npar-1]
                axes[ii,jj]=fig.add_subplot(gs[ii, jj], **axes_kw)
    
    if (_debug_): print(axes)
    
    if (nGroups > 1):
        hist_kw=dict_splicer(hist_kw,nGroups,[1]*nGroups)
        contour_kw=dict_splicer(contour_kw,nGroups,[1]*nGroups)
        scatter_kw=dict_splicer(scatter_kw,nGroups,[1]*nGroups)
        hist2D_kw=dict_splicer(hist2D_kw,nGroups,[1]*nGroups)
    
    
    if (_debug_): print("\nSubplots:")
    for ii in range(npar-1, -1, -1):#(0, npar, 1):
        for jj in range(ii if len(pair_type) == 1 else npar-1, -1, -1):#(0, ii+1, 1):
            pair_num=0 if ii > jj else 1
            
            if (_debug_ == True): print(f" {ii},{jj} ",end='')
            
            if (ii==jj):
                if (nGroups > 1):
                    for kk in range(nGroups):
                        axes[ii,jj].hist(data.iloc[samps][cols[jj][kk]],**hist_kw[kk])
                else:
                    axes[ii,jj].hist(data.iloc[samps][cols[jj]],**hist_kw)
            else:
                if (pair_type[pair_num] == 'scatter'):
                    if (nGroups > 1):
                        for kk in range(nGroups):
                            scatter(data.iloc[samps][cols[jj][kk]],data.iloc[samps][cols[ii][kk]],ax=axes[ii,jj],**scatter_kw[kk])
                    else:
                        scatter(data.iloc[samps][cols[jj]],data.iloc[samps][cols[ii]],ax=axes[ii,jj],**scatter_kw)
                elif (pair_type[pair_num] == 'contour'):
                    if (nGroups > 1):
                        for kk in range(nGroups):
                            contourp(data.iloc[samps][cols[jj][kk]],data.iloc[samps][cols[ii][kk]],ax=axes[ii,jj],**contour_kw[kk])
                    else:
                        contourp(data.iloc[samps][cols[jj]],data.iloc[samps][cols[ii]],ax=axes[ii,jj],**contour_kw)
                elif (pair_type[pair_num] == 'hist2D'):
                    hist2D(data.iloc[samps][cols[jj]],data.iloc[samps][cols[ii]],ax=axes[ii,jj],**hist2D_kw)
                    
            if (labels is not None):
                if (ii == npar-1):
                    axes[ii,jj].set_xlabel(labels[jj])
                if (jj == 0 and ii != 0):
                    axes[ii,jj].set_ylabel(labels[ii])
                if (len(pair_type) == 2 and jj == npar-1 and ii != npar-1):
                    axes[ii,jj].set_ylabel(labels[ii])
                    axes[ii,jj].yaxis.tick_right()
                    axes[ii,jj].yaxis.set_label_position('right')
            
            if (_debug_ == True):
                axes[ii,jj].text(0.1,0.8,f"({ii}, {jj})",transform=axes[ii,jj].transAxes)
            
        if (_debug_ == True): print("")
    
    # turn off redundant tick labeling
    for ii in range(0,npar-1):
        for jj in range(0, ii + 1 if len(pair_type) == 1 else npar):
            axes[ii,jj].xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
            axes[ii,jj].xaxis.offsetText.set_visible(False)
    
    # turn off all but the leftmost columns of each row
    for ii in range(0,npar):
        for jj in range(1, ii):
            axes[ii,jj].yaxis.set_tick_params(which='both',labelbottom=False, labeltop=False)
            axes[ii,jj].yaxis.offsetText.set_visible(False)
    
    # Turn on histogram y-axis labels
    for ii in range(0,npar):
        if (histlabel != None):
            axes[ii,ii].yaxis.tick_right()
            axes[ii,ii].yaxis.set_label_position("right")
            axes[ii,ii].set_ylabel(histlabel)
            axes[ii,ii].tick_params(zorder=3)
        else:
            axes[ii,ii].yaxis.set_tick_params(which='both',labelbottom=False, labeltop=False)
            axes[ii,ii].yaxis.offsetText.set_visible(False)
            
    # Squeeze the axes array
    if (squeeze == True):
        axes=axes.flatten() # flattens the array along along columns then rows
        axes=axes[axes != None] # remove the Nones which fill out the missing corner
        return (fig, axes.item()) if axes.size == 1 else (fig, axes.squeeze())
    else:
        return(fig, axes)


def subplots(naxes=None,nrows=None,ncols=None,va='top',ha='left',wspace=None,hspace=None,
             widths=None,heights=None,sharex='none',sharey='none',squeeze=True,
             figsize=None,axes_kw={},**kwargs):
    """ Adds a set of subplots to figure
    This is a more-generalised wrapper around matplotlib.pyplot.subplot function to allow for irregularly divided grids.
    
    Parameters
    ----------
    naxes : int, optional, default=1
        The number of axes objects to create. The resulting grid formed from specifying naxes is decided by ncols and nrows. The options are:
            * ncols and/or nrows not None:
                Makes sure that naxes can be correctly mapped into the specified grid. If one of nrows or ncols not given, the smallest possible grid will be made.
            * both ncols and nrows are None:
                Decides the best possible grid for this number of axes. Currently, this decision is hard-coded with plans for it to become an automatic decision later.
    
    nrows, ncols : int, optional
        Number of rows/columns of the subplot grid.
    
    va, ha : str, optional, default: 'top', 'left'
        The vertical alignment (va) and horizontal alignment (ha) sets the alignment of grids in the vertical and horizontal directions. Options for ``va`` are 'top', 'bottom' or 'centre' and options for ``ha`` are 'left', 'right' and 'centre'.
    
    wspace, hspace : float, optional
        The horzontal/vertical spacing between figure subplots, expressed as a fraction of the subplot width/height.
    
    width, heights : array-like, optional
        The width/height ratios of the subplot columns/rows expressed as an array of length ncols/nrows.
    
    sharex, sharey : bool or {'none', 'all', 'row', 'col'}, default: False
        As per matplotlib's usage, controls sharing of properties among x (`sharex`) or y (`sharey`)
        axes:
            * True or 'all': x-/y-axis will be shared among all subplots.
            * False or 'none': each subplot x-/y-axis will be independent (default).
            * 'row': each subplot row will share an x-/y-axis
            * 'col': each subplot column will share an x-/y-axis
        When subplots have a shared x-axis along a column, only the x tick labels of the last complete row of the subplot are created. Similarly, when subplots have a shared y-axis along a row, only the y tick labels of the first complete column subplot are created. To later turn other subplots' ticklabels on, use `~matplotlib.axes.Axes.tick_params`.
    
    squeeze : bool, optional, default: True
        As per matplotlib's usage, the following applies:
            * If True, extra dimensions are squeezed out from the returned array of Axes:
                - if only one subplot is constructed (nrows=ncols=1), the
                  resulting single Axes object is returned as a scalar.
                - for Nx1 or 1xM subplots, the returned object is a 1D numpy
                  object array of Axes objects.
                - for NxM, subplots with N>1 and M>1 are returned
                  as a 2D array.
            * If False, no squeezing at all is done: the returned Axes object is always a 2D array containing Axes instances, even if it ends up being 1x1.

        If naxes < ncols*nrows, the only sensible option is to return a 1D numpy array
    
    figsize : 2-tuple of floats, default: rcParams["figure.figsize"] * (ncols, nrows)
        The dimensions of the figure (width, height) in inches. If not specified, the default is to scale
        the default rcParams figure.figsize by the number of rows or columns.

    axes_kw : dict, optional
        Explicit dictionary of kwargs to be parsed to matplotlib `subplot` function.

    **kwargs : Subplot instance properties
        kwargs are used to specify properties of `subplots` instances
        A list of valid `axis` kwargs can be found here:
        [https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes](https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes "Matplotlib.axes.Axes")
    
    Returns
    -------
    fig : pyplot.Figure
    
    axes : pyplot.axes.Axes object or array of Axes objects
        axes may either be a single Axes object or a numpy array if naxes > 1.
        The dimensions of this array are controlled by the squeeze keyword above.
    """
    
    from .base_func import dict_splicer
    
    import warnings
    from matplotlib import rcParams
    from matplotlib.gridspec import GridSpec
    from matplotlib.pyplot import figure, subplot
    from numpy import ceil, array, reshape, empty
    
    gridRef=[[0,0], [1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [2,5], [3,4], [3,4], 
               [4,4], [3,5], [3,5], [4,4], [3,6], [3,6], [4,5], [4,5], [3,7], [5,5], [5,5], [4,6], [5,5]]
    
    # Validate ha and va:
    ha = 'centre' if ha == 'center' else ha
    va = 'centre' if va == 'center' else va
    if ha not in ['left','right','centre']:
        raise ValueError(f"'{ha}' was not recognised for `ha`, must be one of 'left', 'right' or 'centre'")
    if va not in ['top','bottom','centre']:
        raise ValueError(f"'{va}' was not recognised for `va`, must be one of 'top', 'bottom' or 'centre'")

    if (isinstance(sharex, bool)):
        sharex="all" if sharex else "none"
    if (isinstance(sharey, bool)):
        sharey="all" if sharey else "none"
    if (naxes == None): # No number of axes specified
        if (nrows==None): nrows=1
        if (ncols==None): ncols=1
        naxes=nrows*ncols
    else:
        if (ncols == None and nrows == None):
            if (naxes <= 25):
                nrows, ncols=gridRef[naxes] # Get best combination of rows x cols for number of axes
            else:
                raise NotImplementedError(f"The naxes parameter is currently not implemented for naxes > 25.")
        elif (ncols != None and nrows != None):
            if (ncols*nrows != naxes):
                raise ValueError(f"Invalid number of axes ({naxes}) given for number of rows ({nrows}) and columns ({ncols}).") 
        else:
            if (nrows != None):
                ncols=int(ceil(naxes/nrows))
            else:
                nrows=int(ceil(naxes/ncols))
    
    # Number of axes away from filling the gridspec evenly and completely
    delta=(nrows*ncols) - naxes
    
    # Assert that ha/va != center if widths/heights given
    if (widths == None):
        widths=[1]*ncols
    else:
        if ha=='centre': raise ValueError("Cannot set width ratios when ha='centre'")
    
    if (heights == None):
        heights=[1]*nrows
    else:
        if va=='centre': raise ValueError("Cannot set height ratios when va='centre'")
    
    if (figsize==None): # Auto scale default figure size to num cols/rows.
        figsize=(rcParams["figure.figsize"][0]*ncols, rcParams["figure.figsize"][1]*nrows)
    
    fig=figure(figsize=figsize,**kwargs)
    gs=GridSpec(ncols=ncols*2, nrows=nrows*2, hspace=hspace, wspace=wspace,
                  width_ratios=[w for w in widths for kk in range(2)], height_ratios=[h for h in heights for kk in range(2)])
    
    axes_kw=dict_splicer(axes_kw,naxes,[1]*naxes)
    
    # Specify the row/col origin
    row0=(nrows-1)*2 if va=='bottom' else 0
    col0=(ncols-1)*2 if ha=='right' else 0
    axes=empty(naxes, dtype=object) # Create the empty array to hold each axis
    for ii in range(naxes):
        row=2*(nrows-(ii//ncols)-1) if va=='bottom' else 2*(ii//ncols) # current row
        col=2*(ncols-(ii%ncols)-1) if ha=='right' else 2*(ii%ncols) # current column
        
        # Select which axes this axis needs to be shared with
        sharewith={"none": None, "all": axes[0],
                     "row": axes[(ii//ncols)*ncols], "col": axes[ii%ncols]}
        
        axes_kw[ii]["sharex"]=sharewith[sharex]
        axes_kw[ii]["sharey"]=sharewith[sharey]
            
        if (row == (0 if va=='bottom' else (nrows-1)*2) and ha=='centre'):
            axes[ii]=subplot(gs[row:row+2,col+delta:col+delta+2],**axes_kw[ii])
        elif (col == (0 if ha=='right' else (ncols-1)*2) and va=='centre'):
            axes[ii]=subplot(gs[row+delta:row+delta+2,col:col+2],**axes_kw[ii])
        else:
            axes[ii]=subplot(gs[row:row+2,col:col+2],**axes_kw[ii])
    
    # turn off redundant tick labeling
    if sharex in ["col", "all"]:
        if (ha=='centre'):
            warnings.warn("Removing redundant shared xtick labels not possible when ha='centre'")
        else:
            # turn off all but the bottom row
            for ax in axes[ncols:] if va=='bottom' else axes[:naxes-ncols]:
                ax.xaxis.set_tick_params(which='both',
                                         labelbottom=False, labeltop=False)
                ax.xaxis.offsetText.set_visible(False)
    if sharey in ["row", "all"]:
        if (va=='centre'):
            warnings.warn("Removing redundant shared ytick labels not possible when va='centre'")
        else:
            # turn off all but the leftmost column
            for ii, ax in enumerate(axes):
                if (ha=='left'):
                    if (ii%ncols!=0):
                        ax.yaxis.set_tick_params(which='both',
                                                 labelbottom=False, labeltop=False)
                        ax.yaxis.offsetText.set_visible(False)
                elif (ha=='right'):
                    if (ii%ncols+1!=ncols and ii!=naxes-1):
                        ax.yaxis.set_tick_params(which='both',
                                                 labelbottom=False, labeltop=False)
                        ax.yaxis.offsetText.set_visible(False)
    ### Squeeze axes array
    if (squeeze == True):
        # Discarding unneeded dimensions that equal 1.  If we only have one
        # subplot, just return it instead of a 1-element array.
        return (fig, axes.item()) if axes.size == 1 else (fig, axes.squeeze())
    else:
        # Returned axis array will be always 2-d, even if nrows=ncols=1.
        if (naxes != ncols*nrows):
            warnings.warn("squeeze=False not possible when naxes < nrows*ncols.")
            return (fig, axes.squeeze())
        else:
            return (fig, reshape(axes, (nrows, ncols)))
