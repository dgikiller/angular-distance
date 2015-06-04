all:
	mpicc -Wall -lm -ansi "./src/angdist.c" -o angdist
	
clean:
	rm -rf angdist *.o
