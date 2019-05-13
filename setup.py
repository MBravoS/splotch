import atexit
from setuptools import setup
from setuptools.command.install import install

setup(name='splotch',
		version='0.2.7.2',
		description='Simple PLOTs, Contours and Histograms is a small package with wrapper functions designed to simplify plotting calls from matplotlib.',
		url='https://github.com/MBravoS/splot',
		author='Matías A. Bravo Santa Cruz',
		author_email='matias.bravo@icrar.org',
		license='BSD-3-Clause',
		packages=['splotch'],
		#cmdclass={'install':new_install},
		#package_data={'splot/styles':['splot/styles/astro.mplstyle',]},
		install_requires=['numpy>=1.10','matplotlib>=2.2.0','scipy>=1.0.0'],
		zip_safe=False)

