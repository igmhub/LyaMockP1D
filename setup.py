#!/usr/bin/env python
import distutils
from distutils.core import setup

description = "Simple Lya-forest mock maker (independent lines of sight)."

setup(name="lya_mock_p1d", 
      version="0.1.0",
      description=description,
      url="https://github.com/igmhub/LyaMockP1D",
      author="Andreu Font-Ribera",
      py_modules=['lya_mock_p1d'],
      package_dir={'': 'py'})

