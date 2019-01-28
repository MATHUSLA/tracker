# -*- coding: utf-8 -*-
#
# setup.py
#

SETUP_JSON = 'setup.json'

import io
import json
import os
import sys
from shutil import rmtree
from setuptools import setup, find_packages, Command

HERE = os.path.abspath(os.path.dirname(__file__))


def load_json(path, here=HERE):
    """Load JSON from File."""
    with io.open(os.path.join(here, path)) as f:
        return json.load(f)


def get_long_description(path, default='', here=HERE):
    """Get Long Description from README."""
    long_description = default
    with io.open(os.path.join(here, path), encoding='utf-8') as f:
        long_description = '\n' + f.read()
    return long_description


def get_version(path, key='__version__', here=HERE):
    """Get Version from _version File."""
    version = {}
    with open(os.path.join(here, path)) as f:
        exec(f.read(), version)
    return version[key]


class Upload(Command):
    """Upload to PyPI Command."""

    name = 'upload'
    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        """"""

    def finalize_options(self):
        """"""

    def run(self):
        """Run Upload."""
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(HERE, 'dist'))
        except OSError:
            pass
        self.status('Building Source and Wheel (universal) distribution...')
        os.system(f'{sys.executable} setup.py sdist bdist_wheel --universal')
        self.status('Uploading the package to PyPI via Twine...')
        os.system('twine upload --repository-url https://upload.pypi.org/legacy/ dist/*')
        sys.exit()


if __name__ == '__main__':
    about = load_json(SETUP_JSON)

    about['version'] = get_version(about['version_file'])
    del about['version_file']

    about['long_description'] = get_long_description(about['long_description_file'])
    del about['long_description_file']

    if 'exclude' in about['packages']:
        about['packages'] = find_packages(exclude=tuple(about['packages']['exclude']))

    setup(cmdclass={Upload.name: Upload}, **about)
