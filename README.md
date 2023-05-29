# _splotch_: Simple PLOTs, Contours and Histograms

[![Github release](https://img.shields.io/github/release/MBravoS/splotch.svg?label=tag&colorB=54ebff)](https://github.com/MBravoS/splotch/releases) [![PyPI version](https://img.shields.io/pypi/v/splotch.svg?colorB=ff0080)](https://pypi.python.org/pypi/splotch) [![License](https://img.shields.io/pypi/l/splotch.svg)](https://github.com/MBravoS/splotch/blob/master/LICENSE)

<p align="center">
<img src="/example_images/SPLOTCH_logo.png" alt="drawing" width="300"/>
</p>

**_splotch_** is a Python package aimed at simplifying making common plot types used for scientific data visualisation.
The functions and classes included in _splotch_ are a wrappers of varying complexity around several commonly used _matplotlib_ functions and classes.
The core guiding principle in _splotch_ is reducing the number of lines required to produce a given plot, be it by grouping commonly used calls together into a single call or by providing out-of-the-box alternatives to more complex but still common types of plots.
Secondarily, the design and coding style in _splotch_ closely follows _matplotlib_'s conventions, to make the combined use of both packages as pain-free as posible.
Given this guiding principles, _splotch_ is best understood as an add-on to _matplotlib_ and not as a replacement.
_splotch_ also includes the [_SciCM_ colour map package](https://github.com/MBravoS/scicm), from which it sets the Stone colour map as the default colour map.

## Quick start
A simple example of _splotch_ in use:
```python
    import numpy as np, matplotlib.pyplot as plt, splotch as splt
    
    x = np.random.default_rng().normal(size=1000)
    
    splt.hist(x, xlabel='x', ylabel='PDF')
    plt.show()
```

### Example plots made with _splotch_

*TODO: include code snipets for every example.*

 Histograms                 | Scatter Plots
:---:|:---:
| <img src="/example_images/example_hist.png" alt="drawing" width="400"/> |  <img src="/example_images/example_scatter.png" alt="drawing" width="400"/>

 Contour Plots              | Sector Plots             
:---:|:---:
| <img src="/example_images/example_contours.png" alt="drawing" width="400"/> | <img src="/example_images/example_sectorplot.png" alt="drawing" width="400"/>

| Corner Plots              | Subplots                
:---:|:---:
| <img src="/example_images/example_cornerplot.png" alt="drawing" width="400"/>  |  <img src="/example_images/example_subplots.png" alt="drawing" width="400"/>

## Documentation and use guides
_splotch_'s readthedocs page contains the [full documentation](https://splotch.readthedocs.io/en/latest/) of the package.

*TODO: we should consider making an extended quick start guide (as we have in _SciCM_) and/or example pages (as the several example pages in _matplotlib_'s documentation).*

## _splotch_ in the broader colour map Python package ecosystem

*TODO: briefly present how _splotch_ compares to the likes of _plotly_ and _seaborn_*

## Installation guide
The package is available for installation using pip:

    >pip install splotch

Although you may wish to install it directly from GitHub using:

    >pip install git+https://github.com/MBravoS/splotch.git@devel

## How to cite the use of _splotch_
If you are submitting work that uses _splotch_ for publication in a scientific journal, please include a mention of your use.
Some journals include a dedicated section for this purpose (e.g., the [_Software_ section in the Astrophysical Journal](https://journals.aas.org/aastexguide/#software)), which would be the natural place to mention _splotch_ (please include a link to this repository).
If such a section is not included on your journal or choice, please consider adding the following to your acknowledgements:
> The analysis in this work has been performed using the Python programming language, with the open-source package _splotch_ (https://github.com/MBravoS/splotch).

Feel free to expand the previous statement to include the rest of the sofware used in your work!
