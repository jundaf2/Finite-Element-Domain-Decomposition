# Overlapping-Domain-Decomposition-Finite-Element-Method-FE-DDM-of-Electrostatic-Problem
Domain Decomposition Methods have been widely used in the parallel numerical PDE community in the
recent a few decades to tackle large or multi-scale problems with Finite Element Methods. In this course project, three numerical methods with respect to DDM and FEM are carried out to solve one two-dimensional electrostatic problem. The first one FEMCG is using the classic FEM to solve the total solution domain of
the boundary value problem with Conjugated Gradient method for solving the sparse linear system. In the second one FEDDM-ASM, the total solution domain is divided into two-by-two overlapped subdomains and solved using Additive Schwarz Method. In the third one FE-DDM-PCG, the ASM is used as a preconditioner of the Preconditioned Conjugated Gradient method. 

High performance computing technologies provide faster and larger electromagnetic simulations with modern evergrowing HPC hardware. Two of the major goals commercial CAE softwares  pursue are speed-up and
scale-up. The former can be achieved by using multiprocessing, while the latter is by using domain decomposition. The DDM we choose to implement is Additive Schwarz Method based on MPI. The linear system solver we choose is Conjugated Gradient Method using OpenMP. A Matlab version and a serial C++ version are built for verification purpose and a parallel C++ version is implemented to demonstrate the power of parallel numerical algorithms and the compute capability of nowadays parallel machines. Results for both the physical simulation itself and the benchmarking study are shown and analysis in detail.

## Matlab for algorithm verification
The user can directly run the Matlab version with the source code `/matlab project/srcs/Main_Pre_Post_Processing.m`. For the details of the theory related to the math and physics, please refer to `Matlab Verification.pdf`.

## Serial and Parallel C++ versions
The source code of the C++ version is in the `\cplusplus project` folder. Please refer to `C++ Benchmarking.pdf` for more details related to the C++ implementatin and parallel benchmarking.

### Building the cpp source code

To build, create a new directory, anywhere. We will call it "build" (it already exists in this project).
Change to that directory and type

```
cmake <path_to_source_directory>
make
```
in your shell.

This will create programs named
* `build/bin/feddm_serial_cg`
* `build/bin/feddm_parallel_cg`

### Running the cpp programs
Before running, we need to generate the self-defined .txt files that contains the data and information needed by the C++ program using the Matlab scripts 'scripts/autogeogen.m' and 'scripts/autogeogen_forTraditionialSerial.m'. There are already some files generated in the existing 'build' folder.

You can run the serial part by typing
`build/bin/feddm_serial_cg Hmax`

You can run the parallel part by typing
`mpirun -np <nproc> build/bin/feddm_parallel_cg Hmax`

*Hmax* is the positive number that represents the maximum edge length in the mesh, which can reflects the quality and number of unknowns of the mesh. The smaller *Hmax* is, the better the mesh quality is and thus the more unknowns will be solved in the linear matrix equations.

*nproc* is the number of ranks in the MPI comminication procedure and also the number of subdomain. It should be a square equal and greater than 4.