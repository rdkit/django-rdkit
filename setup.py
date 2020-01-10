#!/usr/bin/env python
import os
from setuptools import setup


def get_packages(package):
    """
    Return root package and all sub-packages.
    """
    return [dirpath
            for dirpath, dirnames, filenames in os.walk(package)
            if os.path.exists(os.path.join(dirpath, '__init__.py'))]


setup(
    name='django-rdkit',
    version='0.2.0',
    description='',
    packages = get_packages('django_rdkit'),
    zip_safe=False,
)
