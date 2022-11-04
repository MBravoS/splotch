########################################################################
############################### Color bar ##############################
########################################################################
import warnings
from numpy import arange, argmax, append, array, ceil, empty, full, nanmax, ndarray, reshape, shape
from numpy.random import default_rng as rng
from pandas import DataFrame, Series, RangeIndex
from matplotlib import rcParams
from matplotlib.contour import ContourSet
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import colorbar as mpl_colorbar, figure, gca, sca, subplot
from matplotlib.text import Text
import matplotlib.ticker as tckr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from .base_func import axes_handler, dict_splicer, is_numeric, plot_finalizer

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
    
    # Validate axis input
    old_axis = gca()
    if ax is not None:
        try:  # check if iterable
            _ = (i for i in ax)
            axes = ax
        except (TypeError):
            axes = [ax]
    else:
        axes = [old_axis]
        
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
        raise NotImplementedError("loc must be specified as integer or string. "
                                  "Providing loc as a tuple-like as colorbar anchor position is not yet implemented.")
    
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
        
        cbar = mpl_colorbar(mapper, cax=cax, orientation=orientation, ticks=ticks, **bar_par[ii])
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

    # Return axis to the previous
    sca(old_axis)
    
    return (cbars[0] if len(cbars) == 1 else cbars)
