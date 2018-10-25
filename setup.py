from setuptools import setup

setup(name='splot',
      version='0.1',
      description='SimplePlot is a small package with wrapper functions designed to simplify plotting calls from matplotlib',
      url='https://github.com/MBravoS/splot',
      author='Matías A. Bravo Santa Cruz',
      author_email='matias.bravo@icrar.org',
      license='BSD-3-Clause',
      packages=['splot'],
      install_requires=['numpy','matplotlib','scipy'],
      zip_safe=False)