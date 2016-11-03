Assumes NUMA node 1 is MCDRAM
Requires libnuma and TBB

Usage:

LD_PRELOAD=./libhbm.so HBM_SIZE=<size per process on MB> HBM_THRESHOLD=<threshold in kB> ./binary

Example for 100MB of MCDRAM per process and all allocation over 64K to go to MCDRAM

mpirun -np 72 -env LD_PRELOAD ./libhbm.so -env HBM_SIZE 100 -env HBM_THRESHOLD 64 ./mpihello

