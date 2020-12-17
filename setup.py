import sys
import os
import re
import setuptools
from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "the aattggcc requires setuptools 1.1 higher")

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect version
print(sys.path.insert(0, "scpipelines"))
import version

version = version.__version__

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if major < 3:
    raise SystemExit("""Requires Python 3 or later.""")

aattggcc_packages = find_packages()
aattggcc_package_dirs = {'scpipelines': 'scpipelines'}

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='clusthub',
    version=version,
    description='aattggcc : long-read single-cell nanopore barcode and UMI correction pipeline',
    author='Adam Cribbs',
    author_email='adam.cribbs@ndorms.ox.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='''aattggcc : long-read single-cell nanopore barcode and UMI correction pipeline''',
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="",
    # package contents
    packages=aattggcc_packages,
    package_dir=aattggcc_package_dirs,
    include_package_data=True,
    entry_points={
        "console_scripts": ["aattggcc = scpipelines.entry:main"]
    },
    # other options
    zip_safe=False,
    test_suite="tests",
)
