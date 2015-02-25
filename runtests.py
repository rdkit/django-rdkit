#!/usr/bin/env python

import os, sys

import django
from django.conf import settings
from django.core.management import call_command

DEFAULT_SETTINGS = dict(
    INSTALLED_APPS=(
        #'django_rdkit',
        'tests',
    ),
    DATABASES={
        "default": {
            #'ENGINE': 'django.db.backends.postgresql_psycopg2',
            'ENGINE': 'django_rdkit.backend',
            'NAME': 'test',
            'USER': '',
            'PASSWORD': '',
            'HOST': '',
            'PORT': '',
        }
    },
)


def runtests():
    if not settings.configured:
        settings.configure(**DEFAULT_SETTINGS)

    django.setup()

    parent = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, parent)

    # if the postgresql_psycopg2 backend is used and the 'tests' application
    # is not migrated, it will be synced before the rdkit extension is
    # installed and therefore causing errors.
    #call_command('makemigrations', 'tests')

    #call_command('sqlmigrate', 'tests', '0001')

    sys.exit(call_command('test', 'tests'))

if __name__ == '__main__':
    runtests()
