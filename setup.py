#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='django-rdkit',
      version='0.0.1',
      description='',
      #url='https://github.com/rvianello/django-rdkit',
      packages = find_packages(exclude = [
            'docs',
            'docs.*',
            ]),
     )
