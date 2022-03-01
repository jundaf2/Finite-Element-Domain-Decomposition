#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --job-name proj-benchmark
#SBATCH -p cs
#PBS -S /projects/cs/cs484/sing_shell.sh

export TOTAL_CPUS=$(( SLURM_JOB_NUM_NODES * PBS_NUM_PPN ))

## If not started with PBS, figure out where we are relative to the build directory
##### Snippet from:   http://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
##### End snippet
#IF PBS_O_WORKDIR is not set, we are not running in PBS, choose directory relative to script.
SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}/..}
PBS_O_WORKDIR=${PBS_O_WORKDIR:-${SCRIPT_DIR}/..}

# Moves to the directory the user was in when they ran qsub
cd ${SLURM_SUBMIT_DIR} #assumed to be the source tree

# Check that the script was submit from the right place.
if [ -d "./cmake" ] && [ -d "./writeup" ]
then
    echo "We seem to be in the right place."
else
	echo "Not submit from the right place! Submit from the root of your repository."
	exit 1
fi

# Creates an out-of-tree build directory for CMake and moves to it
mkdir -p ${SLURM_SUBMIT_DIR}/build
pushd ${SLURM_SUBMIT_DIR}/build

:<<EOF
# Bbuild the programs (into the build directory, IE, the current directory)
# then benchmark them. Quit early on failure.
echo "Compiling"
cmake ${SLURM_SUBMIT_DIR} && make
EOF

echo "Testing Serial" >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_serial_more.txt
echo "BEGIN _VARIES_" >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_mpi_openmp_4dm.txt
for FEDDM_RUNSPEC in 1.5 1.4 1.3 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12 0.11 0.1 0.095 0.09 0.085 0.08 0.075 0.07 0.065 
do
	echo "Doing Hmax = ${FEDDM_RUNSPEC} "
	echo "Doing Hmax = ${FEDDM_RUNSPEC} " >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_serial_more.txt
	./bin/feddm_serial_cg ${FEDDM_RUNSPEC} | grep "ASM serial version costs time" >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_serial_more.txt

	#Individual MPI run
	echo "Doing MPI CPUS = 4 "
	echo "Doing MPI CPUS = 4 " >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_mpi_openmp_4dm.txt
	echo "Doing Hmax = ${FEDDM_RUNSPEC} " >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_mpi_openmp_4dm.txt
	srun --mpi=pmi2 --ntasks-per-node ${SLURM_NTASKS_PER_NODE} --cpu-bind=cores --ntasks 4 /projects/cs/cs484/sing_exec.sh ./bin/feddm_parallel_cg ${FEDDM_RUNSPEC} | grep "MPI ASM FEM costs time" >> ${SLURM_SUBMIT_DIR}/writeup/benchmark_mpi_openmp_4dm.txt
done
