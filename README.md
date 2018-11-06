# csolb
Magnetic flux density calculator for solenoid magnets.
<br/>This program can run 100 million calculations of a magnetic field at arbitrary point near a single solenoid in 35 seconds. Stress tests are done with Intel(R) Core(TM) i5-8600 6-core CPU with 16 GB of memory. No GPU is used in the calculation.

### Necessary Libraries
- Intel MKL Library (https://software.intel.com/en-us/mkl)
- OpenMP (Embedded in modern gcc, https://www.openmp.org/)

### Compilation Method
Make sure above necessary libraries are installed in your machine. In GNU programming environment, run make.

### Troubleshooting
- Export your Intel MKL runtime library to LD_LIBRARY_PATH (I provided a bash script of doing it)
- For other issue, please contact <jarin.lee@gmail.com>

### Credit
Programmed by Jaerin Lee <jarin.lee@gmail.com>
<br/>Department of Electrical and Computer Engineering
<br/>Seoul National University

Implementation of core calculation method in the form of matlab code was contributed by Seungyong Hahn, Seoul National University. This work is a redesign of the original code in a more parallel-friendly form, with adding necessary interfaces for ease of use.

### Reference
MW Garrett, "Calculation of Fields, Forces, and Mutual Inductances of Current Systems by Elliptic Integrals," Journal of Applied Physics, vol.34(9), pp. 2567-2573, September 1963.
