#!/bin/bash

# This script loads intel mkl libraries as loadable library.
# You may alter the folder of which the mkl libraries are installed.

export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
