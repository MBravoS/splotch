import atexit
from setuptools import setup
from setuptools.command.install import install

#def _post_install():
#	import splot
#	splot.copy_style()

#class new_install(install):
#	def __init__(self, *args, **kwargs):
#		super(new_install, self).__init__(*args, **kwargs)
#		atexit.register(_post_install)

setup(name='splot',
		version='0.2.2',
		description='SimplePlot is a small package with wrapper functions designed to simplify plotting calls from matplotlib',
		url='https://github.com/MBravoS/splot',
		author='Matías A. Bravo Santa Cruz',
		author_email='matias.bravo@icrar.org',
		license='BSD-3-Clause',
		packages=['splot'],
		#cmdclass={'install':new_install},
		package_data={'splot/styles':['splot/styles/astro.mplstyle',]},
		install_requires=['numpy>=1.10','matplotlib>=2.2.0','scipy>=1.0.0'],
		zip_safe=False)

