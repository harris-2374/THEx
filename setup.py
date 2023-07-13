# -*- coding: utf-8 -*-
import sys
from setuptools import setup, find_packages
from src.thex.version import __version__

if sys.version_info[:2] < (3, 7):
    sys.exit("Tree House Explorer (THEx) requires Python >=3.7 Current Python version: %d.%d" % sys.version_info[:2])

setup(
    name="thex",
    version=f'v{__version__}',
    author="Andrew Harris",
    author_email="ajharris.2374@gmail.com",
    url="https://github.com/harris-2374/THEx",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="Tree House Explorer (THEx) is a novel phylogenomic genome browser.",
    package_dir = {"": "src"},
    packages=find_packages(where="src"),
    package_data={'thex': ['assets/*']},
    scripts = ['src/thex/__main__.py', 'src/thex/app.py', 'src/thexb/__main__.py'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'thex = thex.__main__:main',
            'thexb = thexb.__main__:main',
        ],
    },
)
