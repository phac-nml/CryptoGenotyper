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

exec(open('crypto_typer/version.py').read())

setup(
   name='crypto_typer',
   include_package_data=True,
   version='1.0',
   python_requires='>=3.6',
   setup_requires=['pytest-runner'],
   tests_require=['pytest'],
   packages=find_packages(exclude=['tests','databases']),
   url='https://github.com/christineyanta/crypto_typer',
   author='Christine Yanta',
   author_email='christine.yanta@canada.ca',
   description=('crypto_typer is a tool that allows for subtyping the parasite, Cryptosporidium'),
   keywords='Cryptosporidium typer',
   classfiers=classifiers,
   package_dir={'crypto_typer':'crypto_typer'},

   install_requires=[
      'biopython >=1.70',
      'numpy >=1.1'
   ],

   entry_points={
      'console_scripts': ['crypto_typer=crypto_typer.typer:main'],
   },
)

#crypto_typer -i ./crypto_typer/crypto_typer/example/P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
#python setup.py install --record files.txt && crypto_typer -i crypto_typer/example/P17705_gp60-Crypt14-1F-20170927_gp60F_G07_051.ab1 -m gp60 -t forward -f gp60F -o test
#python setup.py install --record files.txt && crypto_typer -i ./crypto_typer/example/P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
#python setup.py install --record files.txt && crypto_typer -i ./crypto_typer/example/ -m 18S -t contig -f SSUF -r SSUR -o test
#python setup.py install --record files.txt && crypto_typer -i ./test/ -m gp60 -t contig -f gp60F -r gp60R -o test
#crypto_typer -i P17705_Crypto16-2F-20170927_SSUF_G12_084.ab1 -m 18S -t forward -f SSUF -o test
