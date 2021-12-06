"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import re

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='ultrafastFitFunctions',
    version='0.3.0',
    packages=['ultrafastFitFunctions'],
    url='https://github.com/EmCeBeh/ultrafastFitFunctions',  # Optional
    install_requires=['numpy', 'scipy', 'uncertainties'],  # Optional
    license='',
    author='Martin Borchert',
    author_email='martin.b@robothek.de',
    description='Collection of common fit functions for discribing ultrafast and static physics phenomena',  # Required
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)
)