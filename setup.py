# -*- coding: utf-8 -*-
import sys

from setuptools import setup

if sys.version_info[:2] < (3, 7):
    sys.exit("TreeHouseExplorer is only meant for Python >=3.7"
             "Current Python version: %d.%d" % sys.version_info[:2])

setup(
    name="thex",
    version="1.0.0",
    author="Andrew Harris",
    author_email="ajharris.2374@gmail.com",
    url="https://github.com/harris-2374/THEx",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="TreeHouseExplorer is a novel phylogenomic genome browser.",
    package_dir = {"": "src"},
    packages=["thex", "thex.assets", "thex.apps", "thex.apps.utils", "thexb"],
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
