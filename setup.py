#!/usr/bin/env python

from distutils.core import setup
setup(name='spectra',
      version='0.2',
      description='Basic spectral analysis using lmfit',
      author='Rohan Isaac',
      author_email='rohan_isaac@yahoo.com',
      url='https://github.com/rohanisaac/spectra',
      packages=['spectra'],
      install_requires=['numpy', 'lmfit', 'scipy'],
      )
