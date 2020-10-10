#!/usr/bind/env python3
import os
from distutils.core import setup
from setuptools import find_packages

author='Christine Yanta'

classifiers = """
Development Status :: Version 1.0
Intended Audience :: Science/Research
Topic :: Science :: Bioinformatics
Programming Language :: Python :: 3.6
""".strip().split('\n')

def read(fname):
   return open(os.path.join(os.path.dirname(__file__), fname)).read()

exec(open('CryptoGenotyper/version.py').read())

setup(
   name='cryptogenotyper',
   include_package_data=True,
   version='1.0',
   python_requires='>=3.6',
   setup_requires=['pytest-runner'],
   tests_require=['pytest'],
   packages=find_packages(exclude=['tests']),
   url='https://github.com/phac-nml/CryptoGenotyper',
   author='Christine Yanta',
   author_email='christine.yanta@canada.ca',
   description=('CryptoGenotyper is a tool that allows for subtyping the parasite, Cryptosporidium'),
   keywords='Cryptosporidium typer',
   classfiers=classifiers,
   package_dir={'CryptoGenotyper':'CryptoGenotyper'},

   install_requires=[
      'biopython >=1.70,<1.78',
      'numpy >=1.1'
   ],

   entry_points={
      'console_scripts': ['cryptogenotyper=CryptoGenotyper.typer:main'],
   },
)
