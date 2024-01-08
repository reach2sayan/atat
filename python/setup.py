from distutils.core import setup, Extension
from setuptools import find_packages

setup(
    name='ATATPythonFunctions',
    version='0.1.0',
    author='Sayan Samanta',
    author_email='sayan_samanta@brown.edu',
    packages=find_packages(),
    scripts=['../python/pythonfunctions.py'],
    install_requires=[
        "ase >= 3.22.1",
        "numpy >= 1.23.5",
        "scipy >= 1.9.0",
        "sympy >= 1.11",
        "matplotlib >= 3.5.3",
        "pandas >= 1.4.3",
    ],
)
