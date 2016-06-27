#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "computeFlux.h"
#include "macroDefinition.h"
#include "updateFlux.h"
#include "visual.h"
#include "computeTimeStep.h"

# define physx( i, mex, NX ) ( (i) + mex*NX )	// @ANDY:added
# define physy( j, mey, NY ) ( (j) + mey*NY ) 	// @ANDY:added

void printMatrixf(int n_grid, double *matrix)
{
	for (int x = 0; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			printf("matrix[%d][%d]: %lf\t",x,y, matrix[x*n_grid + y] );
		}
		printf("\n");
	}
}

int main(int argc, char **argv)
{
    /* initialise MPI */
    int x, y;				// @ANDY:added
    int rank, size;			// @ANDY:added
	int mex, mey;  // (0,0) in (x,y), is at the bottom left of the grid of tiles, whre tiles are processors
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	/* initial value set up */
	int cellsize = 1;
	int n_grid = 50;
	int length = n_grid * cellsize;
	int totalNumberofTimeStep = 1000;
	int plottingStep = 1;
	int LTS_levels = 5;  /* LTS level: GTS: 1, LTS: 2,3,... sets the numbr of local dt values (resolution) applied to the grid */
	int level_A, level_B; 			// @ANDY:LTS
	int idum;  						// @ANDY:LTS
    printf("IM rank: %d", rank);
	double dt0 = 0.1; 				// @ANDY:LTS: @dt0 is just default, cahnged l8er?
	// double dt_dx = dt0/cellsize; // @ANDY:LTS
	double crmax; 					// @ANDY:LTS crmax is returned by computeTimeStep.c

	const char *szProblem;
	szProblem = "result";

	/* Initialisation & memory allocation */
	// double amax;
	double *h, *u, *v, *F, *G, *U, *dt, *lambdaf, *lambdag;
	int *levelf, *levelg, *levelc; // levelc is most important: its elements are u/ to determine whether some grid point needs to be updated at any given time (iteration). Grid points that are distant to the shockwave experience less change in height/velocity/etc, so those will be updated less frequently (larger dt) than other grid points.

	h = malloc(n_grid*n_grid*sizeof(double));
	u = malloc(n_grid*n_grid*sizeof(double));
	v = malloc(n_grid*n_grid*sizeof(double));
	dt = malloc(n_grid*n_grid*sizeof(double));
	lambdag = malloc((n_grid + 1)*n_grid*sizeof(double)); 	// @ANDY:LTS: arrays compared to GTS
	lambdaf = malloc((n_grid + 1)*n_grid*sizeof(double)); 	// @ANDY:LTS
	F = malloc((n_grid+1)*n_grid*3*sizeof(double));
	G = malloc((n_grid+1)*n_grid*3*sizeof(double));
	U = malloc(n_grid*n_grid*3*sizeof(double));

	levelc = malloc(n_grid*n_grid*sizeof(int)); 			// @ANDY:LTS
	levelf = malloc((n_grid+1)*n_grid*sizeof(int)); 		// @ANDY:LTS
	levelg = malloc((n_grid+1)*n_grid*sizeof(int));			// @ANDY:LTS

	/* initialisation */
	for (int x = 0; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			h[x*n_grid + y] = 0.1;
			u[x*n_grid + y] = 0.0;
			v[x*n_grid + y] = 0.0;

			U[ (x*n_grid + y)*3] = h[x*n_grid + y];
			U[ (x*n_grid + y)*3 + 1] = u[x*n_grid + y] * h[x*n_grid + y];
			U[ (x*n_grid + y)*3 + 2] = v[x*n_grid + y] * h[x*n_grid + y];
			levelc[x*n_grid + y] = LTS_levels*1;
			dt[x*n_grid + y] = dt0;  // @dt0: @Q:xiao: why is dt0 always the same? I cannot see where in the code it changes to allow multiple timescale?
		}
	}
	/*initialise h*/
	for (int x = 0; x < 6; ++x)
	{
		for (int y = n_grid - 6; y < n_grid; ++y)
		{
			h[x*n_grid + y] = 1.0;
			U[ (x*n_grid + y)*3] = h[x*n_grid + y];
		}
	}

	for (int x = 0; x < n_grid + 1; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			levelf[x*n_grid + y] = LTS_levels*1;
			lambdaf[x*n_grid + y] = 0.0;
			for (int i = 0; i < 3; ++i)
			{
				F[ (x*n_grid + y)*3 + i ] = 0.0;
			}
		}
	}
	for (int x = 0; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid + 1; ++y)
		{
			levelg[x*(n_grid+1) + y] = LTS_levels*1;
			lambdag[x*(n_grid+1) + y] = 0.0;
			for (int i = 0; i < 3; ++i)
			{
				G[ (x*(n_grid+1) + y)*3 + i ] = 0.0;
			}
		}
	}

	/* start simulation */
	level_B = LTS_levels;
	for (int i = 1; i <= totalNumberofTimeStep; ++i)
	{

		level_A = level_B;
		/* initialise different level */
		if ((i-1)%2 == 0)
		{
			level_B = 1;
		}else{
			for (int j = 1; j <= LTS_levels; ++j)
			{
				idum = (int)pow(2,j-1);
				if (i%idum == 0)
				{
					level_B = j;
				}
			}
		}
		// printf("level_A: %d\n", level_A);
		/* compute fluxes*/
		computeFlux(U, F, G, levelf, levelg, lambdaf, lambdag, n_grid, level_A);

		/* Assign LTS levels  & calculate timestepping*/
		if (level_A == LTS_levels)
		{
			calculateTimeStep(lambdaf, lambdag, dt, levelc, levelf, levelg, LTS_levels, cellsize, 
				n_grid, dt0, &crmax);
		}
		
		// for (int x = 0; x < n_grid + 1; ++x)
		// {
		// for (int y = 0; y < n_grid; ++y)
		// {
		// 	printf("levelf[%d][%d]: %d\t",x,y, levelf[x*n_grid + y] );
		// }
		// printf("\n");
		// }
		// for (int x = 0; x < n_grid; ++x)
		// {
		// for (int y = 0; y < n_grid + 1; ++y)
		// {
		// 	printf("levelg[%d][%d]: %d\t",x,y, levelg[x*(n_grid+1) + y] );
		// }
		// printf("\n");
		// }
		// for (int x = 0; x < n_grid; ++x)
		// {
		// for (int y = 0; y < n_grid; ++y)
		// {
		// 	printf("levelc[%d][%d]: %d\t",x,y, levelc[x*n_grid + y] );
		// }
		// printf("\n");
		// }
		// for (int x = 0; x < n_grid; ++x)
		// {
		// for (int y = 0; y < n_grid; ++y)
		// {
		// 	printf("dt[%d][%d]: %lf\t",x,y, dt[x*n_grid + y] );
		// }
		// printf("\n");
		// }
		// printMatrixf(n_grid, levelf);


		/* updating the fluxes*/
		updateFlux(U, F, G, levelc, n_grid, dt, level_B);
	// for (int i = 0; i < 3; ++i)
	// {
	// 	for (int x = 0; x < n_grid; ++x)
	// 	{
	// 		for (int y = 0; y < n_grid; ++y)
	// 		{
	// 			printf("U: %lf \t", U[(x*n_grid + y)*3 + i]);
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }

		
		//printf("Time Step = %d, amax = %lf \n", i, amax);
		printf("Time Step = %d, Courant Number = %lf, level_B: %d \n", i, crmax, level_B);
		// /* write vtk file*/
		// if(i % plottingStep == 0){
  //       	write_vtkFile(szProblem, i, length, n_grid, n_grid, cellsize, cellsize, U, dt);
  //       }
	}

	/* memory deallocation */
	free(h);
	free(u);
	free(v);
	free(F);
	free(G);
	free(U);
	free(levelf);
	free(levelc);
	free(levelg);
	free(lambdag);
	free(lambdaf);
    MPI_Finalize();
}
