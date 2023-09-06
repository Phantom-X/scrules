# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/4 7:46
@Email: 2909981736@qq.com
"""

from setuptools import setup, find_packages
__version__ = '1.0.2'

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name="scrules",
    version=__version__,
    author="phantom",
    author_email="2909981736@qq.com",
    description="A library that mines association rules between genes and mines regulatory networks through the association rule algorithm fp-growth",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/PH-Open/scrules",
    packages=find_packages(),
    install_requires=[
        "pyfim",
        "igraph",
        "numpy",
        "pandas",
        "matplotlib",
        "pymysql",
        "sqlalchemy"
    ],
    classifiers=[
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        'Programming Language :: Python :: Implementation :: CPython',
    ],
)
