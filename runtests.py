#!/usr/bin/env python

import os, sys

import django
from django.conf import settings
from django.core.management import call_command

DEFAULT_SETTINGS = dict(
    INSTALLED_APPS=(
        'django_rdkit',
        'tests',
    ),
    DATABASES={
        "default": {
            'ENGINE': 'django.db.backends.postgresql_psycopg2',
            'NAME': os.environ.get('DJANGO_DB', 'test'),
            'USER': os.environ.get('DJANGO_USER', ''),
            'PASSWORD': os.environ.get('DJANGO_PASSWORD', ''),
            'HOST': os.environ.get('DATABASE_HOST', ''),
            'PORT': os.environ.get('DATABASE_PORT', ''),
        }
    },
    DEFAULT_AUTO_FIELD='django.db.models.BigAutoField',

    DJANGO_RDKIT_MOL_SERIALIZATION=os.environ.get('DJANGO_RDKIT_MOL_SERIALIZATION', 'BINARY'),
)


def runtests():
    if not settings.configured:
        settings.configure(**DEFAULT_SETTINGS)

    django.setup()

    parent = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, parent)

    sys.exit(call_command('test', 'tests', *sys.argv[1:]))

if __name__ == '__main__':
    runtests()
