Description for Xnet:

To perform the nuclear reaction calculations in our D6 models, we will make use of the Xnet integrated FLASH code. Xnet is a nuclear reaction code that evolves the abundances of nuclides based on the reaction rates. It is a modular code, designed for arbitrary size networks of astrophysical interests. It makes uses of Bader-Deuflhard method to solve a system of first-order differential equations. It is written in FORTRAN95 and makes uses of MPI, OpenACC, and OpenMP to parallelize across multiple GPUS and cores across numerous nodes.

These are all the things relavant to xnet:

- The underlying mathematical formulation (e.g., ODE, PDE).
  - First Order Differential Equation.

- Particular libraries required by the simulation and analysis software, algorithms and numerical techniques employed (e.g., finite element, iterative solver), programming languages, and other software used.
  - Libraries:
    - LAPACK
    
  - numerical techniques employed:
    - BDF integrator
    - Jacobi iterative method
  - Programming languages:
    - Fortran
    - Python
    - C/C++
  - other software:
  
- Parallel programming model(s) used (e.g., MPI, OpenMP, Pthreads, CUDA, OpenACC).
  - MPI, OpenACC, CUDA, OpenMP
 
- Project workflow, including the role of analysis and visualization and the availability of checkpoint and restart files.

- I/O requirements (e.g., amount, size, bandwidth, etc) for restart, analysis, and workflow. Highlight any exceptional I/O needs.
