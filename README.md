# SingleFileHardDisks - Structural Properties

This is the source code accompanying [1] to compute the structural properties (nth-neighbor probability distribution function, radial distribution function and static structure factor) of the single-file hard disks fluid, as described in [1].

## :ferris_wheel: Dependencies

The source code is provided as a C++ project and has the following prerequisites
- gcc (tested in gcc 11.2.0) or a compatible compiler.
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library
- Make
- git

## :classical_building: Compilation

The project is provided ready to be compiled under Linux platforms
1. `` git clone  https://github.com/amonterouex/SingleFileHardDisks.git``
2. `` cd SingleFileHardDisks ``
4. Make sure the location of your local Eigen library matches the standard one provided in line
`` EIGEN = /usr/include/eigen3/``
of the Makefile. If not, modify that line to point at the correct Eigen location.
3. `` make ``


## :airplane: Execution

Two different files are needed to correctly run the code

### input.dat

Contains information about the system and the structural properties one wishes to compute. User can modify the values in this file (that will be read later on by the main program).

### main

The correct compilation of the source code should produce an executable named main. The programs reads file "input.dat", computes the required quantities, and writes the output to a file.


## :envelope_with_arrow: Contact

If you have trouble compiling or using this software, if you found a bug or if you have an important feature request, you may contact us at <anamontero@unex.es>

## :books: References
[1] Ana M. Montero and A. Santos, Structural Properties of Hard-Disk Fluids under Single-File Confinement, 	[arXiv:2304.14290](https://arxiv.org/pdf/2304.14290)
