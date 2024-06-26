# -*- coding: utf-8 -*-
"""
Setup file for FLOW_API
"""
import os
from setuptools import setup, find_packages

setup(name='flow_api',
      version='v0.1',
      description='FLOW API',
      url='https://gitlab.windenergy.dtu.dk/eu-flow/wp4/FLOW_API/',
      author='EU FLOW consortium',
      author_email='juqu@dtu.dk',
      packages=find_packages(),
      install_requires=[
          'pytest',  # for testing
          'pytest-cov',  # for calculating coverage
          'pycodestyle',  # Check code style
          'py_wake'
      ],
      zip_safe=True)
