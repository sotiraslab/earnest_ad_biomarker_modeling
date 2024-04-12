#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:22:42 2024

@author: tom.earnest
"""

from os import path

from setuptools import setup

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
    
requirements = [
    'matplotlib',
    'numpy',
    'pandas',
    'pingouin',
    'scipy',
    'scikit-learn',
    'statsmodels'
    ]

setup(name='atn_modeling',
      version='0.0.1',
      description='ATN modeling modules.',
      author='Tom Earnest',
      author_email='tom.earnest@wustl.edu',
      license='MIT',
      packages=['atn_modeling'],
      install_requires=requirements,
      long_description=long_description,
      long_description_content_type='text/markdown')