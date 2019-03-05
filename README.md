# MP-CBM-Z-v1.0
MP CBM-Z v1.0:
This tar bag includes the baseline and optimized CBM-Z codes and a 3D test case with 160*148*20 grid boxes.
The main directory contains:
data/                  Directory for the test case data
cbmz_base/             Directory for the baseline CBM-Z codes
cbmz_opt/              Directory for the optimized CBM-Z (MP CBM-Z) codes.
cbmz_mpi/              Directory for the MPI version MP CBM-Z codes
cbmz_openmp/           Directory for the OpenMP version MP CBM-Z codes
bin/                   Directory for executable files

Compile and test the codes with Intel Compiler:
1. Getting into the corresponding code directory(base, opt, openmp and mpi), and using `make` to generate the binary file;
   e.g.  cd ./cbmz_base
         make
   After compiling codes, the executable files were put into the bin direcotry.

2. Runing the binary file in the bin factory, and the test case would be read automatically. The test would be repeated 10 times and the consuming times would be prtinted after running:
   e.g,  ./cbmz_opt.cpu.fast
  
3. Runing the MPI version MP CBM-Z, please use the command as follow:
   e.g. mpirun -n 4 ./cbmz_opt.mpi
   This command would use 4 MPI processes to run the MP CBM-Z.
4. Runing the OpenMP version MP CBM-Z, please use the command as follow:
   e.g. export OMP_NUM_THREADS=4
        export OMP_STACKSIZE=32M
	./cbmz_opt.openmp
   These command would use 4 OpenMP threads to run the MP CBM-Z.
