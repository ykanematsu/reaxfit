#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import find_packages, setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
        name="reaxfit",
        author="ykanematsu",
        use_scm_version=True,
        setup_requires=["setuptools_scm"],
        license='MIT',
        url='https://github.com/ykanematsu/reaxfit',
        description='parameter fitting for ReaxFF',
        long_description=long_description,
        long_description_content_type='text/markdown',
        #install_requires=["scipy","numpy","lammps"],
        package_dir={"": "src"},
        packages=find_packages(where="src"),
)        
