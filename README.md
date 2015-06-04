# Angular Distance
---
Angular Distance is program for calculating distances between a given set of stars.<br>
This program uses MPI for calculating stars' distance with multiple process, fairly splitting work between these.
## Compilation & Running
Compile the source with:
```
mpicc -Wall -lm -ansi angdist.c -o angdist
```
or using the Makefile.<br>
After this you can run it with:
```
mpirun -n <processes> ./angdist <K-intervals> <stars-number> <input-file> [-debug|-gnuplot]
```
## License 
GPLv3
## Author
Francesco Ferro