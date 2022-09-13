import atexit
from setuptools import setup
from setuptools.command.install import install

setup(name='splotch',
      version='0.5.10.1',
      description='Simple PLOTs, Contours and Histograms is a small package with wrapper functions designed to simplify plotting calls from matplotlib.',
      url='https://github.com/MBravoS/splot',
      author='Matías A. Bravo Santa Cruz',
      author_email='matias.bravo@icrar.org',
      license='BSD-3-Clause',
      packages=['splotch'],
      package_dir={'splotch': 'src/splotch'},
      include_package_data=True,
      package_data={'src':['styles/*.style']},
      #cmdclass={'install':new_install},
      install_requires=['numpy>=1.15','matplotlib>=3.0.0','scipy>=1.1.0','sympy>=1.2'],
      extras_require={ # optional packages
                      'astro':['astropy']
      },
      zip_safe=False)

