This is a repository for the MPI-based C++ code to generate a phase diagram for the Domany-Kinzel-like model of a two-dimensional range expansion. 


Source code file list

- `phase.cpp`: main program
- `rancombi.cpp`: template definition for random number generators 
- `randomc.h`: header file for random number generators
- `ranrotw.cpp`: random number generator RANROT type W implementation
- `ranrotb.cpp`: random number generator RANROT type B implementation


The code is complied via "mpic++ phase.cpp".  The random number generator files must be included for proper linking.  The local configuration must also support the MPI libraries. Please read the "DOCUMENTATION" file for details of implementation.
