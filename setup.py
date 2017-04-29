#! /usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="kdtable",
    version='v0.1',
    author="Johnny Greco",
    author_email="jgreco@astro.princeton.edu",
    packages=["kdtable"],
    url="https://github.com/johnnygreco/Johnny Greco",
    license="MIT",
    description="astropy table with kdtree property for querying large catalogs",
)
