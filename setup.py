#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import find_packages, setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
        name="reaxfit",
        version="0.0.1",
        install_requires=["scipy","numpy","lammps"],
        package_dir={"": "src"},
        packages=find_packages(where="src"),
)        
