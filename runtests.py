#!/usr/bin/env python

import os, sys

import django
from django.conf import settings
from django.test.runner import DiscoverRunner


DEFAULT_SETTINGS = dict(
    INSTALLED_APPS=(
        'django_rdkit',
        'tests',
    ),
    DATABASES={
        "default": {
            'ENGINE': 'django_rdkit.db.backends.postgresql_psycopg2',
            'NAME': 'test',
            'USER': '',
            'PASSWORD': '',
            'HOST': 'localhost',
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

    failures = DiscoverRunner(
        verbosity=1, 
        interactive=True, 
        failfast=False
    ).run_tests(['tests'])

    sys.exit(failures)


if __name__ == '__main__':
    runtests()
