#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "attrs",
    "requests",
    "cachetools",
]

test_requirements = [
    "pytest>=3",
]

setup(
    author="David James Skelton",
    author_email="james.skelton@newcastle.ac.uk",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="A Python library for interfacing with the PNeDRex API",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="python_nedrex",
    name="python_nedrex",
    packages=find_packages(include=["python_nedrex", "python_nedrex.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/james-skelton/python_nedrex",
    version="0.1.1",
    zip_safe=False,
)
