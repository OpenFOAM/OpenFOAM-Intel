## Scalable motorbike benchmark

There are scripts to configure a motorbike example with minimal effort.

First decide the size and mesh.  This is configured with the inital blockmesh size:

    ./Mesh <X> <Y> <Z>

Here is an example of the sizes of mesh that are created for different values:

|   X   |   Y   |   Z   |  MCells  |
|-------|-------|-------|----------|
|    20 |     8 |     8 |     0.35 |
|    60 |    24 |    24 |     5.38 |
|    80 |    32 |    32 |    11.2  |
|    90 |    36 |    36 |    15.5  |
|   100 |    40 |    40 |    20    |

It is recommended to keep the same proportions here.

Note: All meshing will be run with 16 procs although you can update if necessary.

Now, you can set up the case for a given number of processes (calling setup again will wipe existing partitions):

    ./Setup <NPROCS>

This will decompose and run potentialFoam.

Now, all you need to do is solve.  Either run simpleFoam or call:

    ./Solve


