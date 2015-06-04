/* angdist.c
 *
 * A simple MPI program that calculate all angular distance between stars of a given set.
 *
 * Compile command:
 * mpicc -Wall -lm -ansi angdist.c -o angdist
 *
 * Run command:
 * mpirun -n <processes> ./angdist <K intervals> <stars number> <input file> [-debug|-gnuplot]
 *
 * Author: Francesco Ferro <skullbocks@dark-lab.net>
 * License: GPLv3
 *
 */

/********************************************************************/
/*                           Library Import                         */
/********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

/********************************************************************/
/*                    	  Costants definition                       */
/********************************************************************/

#define DIMENSIONS 3
#define MASTER 0
#define N_ELEM_STAR 1
#ifndef M_PI
#define M_PI           3.14159265358979323846264338327
#endif

/********************************************************************/
/*                   	 Structures Definition                      */
/********************************************************************/

/* Star structure */
typedef struct star
{
	double p[DIMENSIONS];
} star_t;

/* Timer structure */
typedef struct {
	double tstart;
	double tstop;
} mytimer_t;

/********************************************************************/
/*                   	Functions definition                        */
/********************************************************************/

/* Start timer function */
void timer_start( mytimer_t* t )
{
	t->tstart = MPI_Wtime();
}

/* Stop timer function */
void timer_stop( mytimer_t* t )
{
	t->tstop = MPI_Wtime();
}

/* Time elapsed function */
double timer_elapsed( const mytimer_t* t )
{
	return t->tstop - t->tstart;
}

/* Stars acquisition from file */
short int acquire_stars(char *file_name,
                        unsigned int stars_number,
                        star_t *stars,
                        unsigned short debug,
                        unsigned short gnuplot)
{
	unsigned int i,
	         j;
	short int status = 0;
	FILE *stars_file;
	stars_file = fopen(file_name, "r");
	if (stars_file == NULL)
	{
		printf("Unable to open the file selected %s\n",
		       file_name);
		status = 1;
	} else {
		for (i = 0; i < stars_number; i++)
		{
			for (j = 0 ; j < DIMENSIONS; j++)
			{
				fscanf(stars_file, "%lf", &stars[i].p[j]);
				if (debug == 1 && gnuplot == 0) {
					if (j == 0)
						printf("%u\t", i);
					printf("%lf\t", stars[i].p[j]);
				}
			}
			if (debug == 1 && gnuplot == 0)
				printf("\n");
		}
		fclose(stars_file);
	}
	return status;
}

/* Calculate row index from array index */
unsigned int row_index( unsigned long array_index,
                        unsigned int matrix_size )
{
	long double m = matrix_size - 1;
	long double r_index = (-2 * m - 1 + sqrt( (4 * m *
	                       (m + 1) - 8 * (long double)array_index - 7) ))
	                      / -2;
	if ( r_index == (long double)(unsigned int) r_index )
		r_index--;
	return (unsigned int) r_index;
}

/* Calculate colum index from array index and row index */
unsigned int column_index( unsigned long array_index,
                           unsigned int r_index,
                           unsigned int matrix_size )
{
	return (unsigned int) (array_index -
	                       ((unsigned long)matrix_size - 1) *
	                       (unsigned long)r_index +
	                       (unsigned long)r_index *
	                       ((unsigned long)r_index + 1) / 2) + 1;
}

/* Angular distance function */
double compute_ang_dist(star_t star_i,
                        star_t star_j)
{
	unsigned short i;
	double numerator = 0.0,
	       coeff_i_sum = 0.0,
	       coeff_j_sum = 0.0;
	for (i = 0; i < DIMENSIONS; ++i)
	{
		numerator += star_i.p[i] * star_j.p[i];
		coeff_i_sum += star_i.p[i] * star_i.p[i];
		coeff_j_sum += star_j.p[i] * star_j.p[i];
	}
	return acos(numerator / sqrt(coeff_i_sum * coeff_j_sum));
}
/* Main Function */
int main(int argc, char *argv[])
{
	int tasks_number,
	    taskid,
	    status = 0;
	double intervals_size,
	       tmp_distance = 0.0;
	unsigned long i,
	         start_distance = 0,
	         end_distance = 0,
	         distance_number = 0, /*N(N - 1)/2*/
	         distance_number_local = 0,
	         distance_number_prev_p = 0;
	unsigned int *k_intervals,
	         *k_intervals_local,
	         stars_number = 0,
	         k_intervals_number = 0,
	         j,
	         star_i = 0,
	         star_j = 0,
	         tmp_interval = 0,
	         gnuplot = 0,
	         debug = 0;
	long tmp_star_number = 0,
	     tmp_k_number = 0;
	star_t *stars = NULL;
	mytimer_t t;
	/* MPI inizialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &tasks_number);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	/* MPI structure for stars inizialization */
	int l_elem_star[N_ELEM_STAR] = {DIMENSIONS};
	MPI_Datatype tipi[N_ELEM_STAR] = {MPI_DOUBLE};
	MPI_Datatype mpi_star;
	MPI_Aint offsets[N_ELEM_STAR];
	offsets[0] =  offsetof(star_t, p);

	MPI_Type_create_struct(N_ELEM_STAR,
	                       l_elem_star,
	                       offsets,
	                       tipi,
	                       &mpi_star);
	MPI_Type_commit(&mpi_star);

	if (argc == 4 || argc == 5) {

		tmp_star_number = atol(argv[2]);
		tmp_k_number = atol(argv[1]);
		if (argc == 5)
		{
			/* Checking extra flags */
			if (strcmp(argv[4], "-gnuplot") == 0)
			{
				gnuplot = 1;
			}
			if (strcmp(argv[4], "-debug") == 0)
			{
				debug = 1;
			}
		}

		/* Abort execution if K or stars number is <= 0 */
		if (tmp_star_number <= 0 || tmp_k_number <= 0)
		{
			if (taskid == 0) {
				printf("Aborting execution, k intervals or stars number is <= 0\n");
				MPI_Abort( MPI_COMM_WORLD, -1 );
			}
		} else {
			stars_number = (unsigned int) tmp_star_number;
			k_intervals_number = (unsigned int) tmp_k_number;
		}
		/* Total number of distance we need to compute */
		distance_number = ((unsigned long)stars_number *
		                   ((unsigned long)stars_number - 1)) / 2;
		stars = (star_t *) calloc(stars_number, sizeof(star_t));
		if (taskid == MASTER) {
			/* Printing program heading */
			if (gnuplot == 0)
			{
				printf("#########################################################\n");
				printf("##                                                     ##\n");
				printf("##                                                     ##\n");
				printf("##                                                     ##\n");
				printf("##      This program calculate the distribution of     ##\n");
				printf("##           angular distance between stars            ##\n");
				printf("##                                                     ##\n");
				printf("##                        by                           ##\n");
				printf("##                  Francesco Ferro                    ##\n");
				printf("##                                                     ##\n");
				printf("##                                                     ##\n");
				printf("#########################################################\n");
			}
			/* Stars acquisition */
			status = acquire_stars(argv[3], stars_number, stars, debug, gnuplot);
			if (status == 0) {
				if (debug == 1 && gnuplot == 0) {
					printf("Total distance number: %lu\n", distance_number);
					printf("Total star number: %u\n", stars_number);
				}
				k_intervals = calloc(k_intervals_number, sizeof(unsigned int));
			} else {
				printf("Error reading stars file\n");
				MPI_Abort( MPI_COMM_WORLD, -1 );
			}
		}
		k_intervals_local = calloc(k_intervals_number, sizeof(unsigned int));
		distance_number_local = distance_number / tasks_number +
		                        (taskid < distance_number % tasks_number );
		if (taskid > 0) {
			distance_number_prev_p = distance_number / tasks_number +
			                         (taskid - 1 < distance_number
			                          % tasks_number );
			/* Initial distance array index for this process */
			start_distance = (taskid - 1) * distance_number_prev_p +
			                 distance_number_prev_p;
		}
		/* Last distance array index for this process */
		end_distance = start_distance + distance_number_local;
		intervals_size = M_PI / (double) k_intervals_number;
		/* Starting timer */
		timer_start( &t );
		/* Stars broadcast to all process */
		MPI_Bcast(stars,            /* Message */
		          stars_number,  /* Dimension */
		          mpi_star,        /* Data Type */
		          MASTER,            /* Sender */
		          MPI_COMM_WORLD);
		/* Computing the angular distance distribution */
		for (i = start_distance; i < end_distance; ++i)
		{
			star_i = row_index(i,
			                   stars_number);
			star_j = column_index(i,
			                      star_i,
			                      stars_number);
			tmp_distance = compute_ang_dist(stars[star_i],
			                                stars[star_j]);
			tmp_interval = (unsigned int) (tmp_distance / intervals_size);
			k_intervals_local[tmp_interval]++;
		}
		/* Reducing the distribution array */
		MPI_Reduce(k_intervals_local,      /* Send buffer */
		           k_intervals,             /* Reception buffer */
		           k_intervals_number,      /* Number of buffer elements */
		           MPI_INT,      			 /* Data type */
		           MPI_SUM,         		 /* Reduction operation */
		           MASTER,          		 /* Destination process */
		           MPI_COMM_WORLD);
		/* Timer stop */
		timer_stop( &t );
		if (taskid == MASTER) {
			if (gnuplot == 0)
			{
				printf("K\tNumber of elements\n");
			}
			for (j = 0; j < k_intervals_number; ++j)
			{
				printf("%lf\t%u\n", (j + 1)*intervals_size,
				       k_intervals[j]);
			}
			free(k_intervals);
			printf("K interval size: %lf \n",
			       intervals_size);
			printf("Total distance calculated: %lu \n",
			       distance_number);
			printf("Time elapsed for parallel operations: %lf s\n",
			       timer_elapsed( &t ) );
		}
		if (debug == 1 && gnuplot == 0) {
			MPI_Barrier(MPI_COMM_WORLD);
			printf("id %d start distance %lu end distance %lu ",
			       taskid,
			       start_distance,
			       end_distance);
			printf("Number of local distances %lu\n", distance_number_local);
			printf("Time elapsed for the process %d: %lf s\n",
			       taskid,
			       timer_elapsed( &t ) );
		}
		free(stars);
		free(k_intervals_local);
	} else {
		if (taskid == MASTER)
		{
			printf("Wrong number of arguments, try:\n");
			printf("mpirun -n <processes> ./angdist <K intervals> <stars number> <input file>\n\n");
			printf("For more info use:\n");
			printf("mpirun -n <processes> ./angdist <K intervals> <stars number> <input file> -debug\n");
			printf("For output gnuplot friendly:\n");
			printf("mpirun -n <processes> ./angdist <K intervals> <stars number> <input file> -gnuplot\n");
			MPI_Abort( MPI_COMM_WORLD, -1 );
		}
	}
	/* MPI finalization */
	MPI_Finalize();
	return (status);
}
