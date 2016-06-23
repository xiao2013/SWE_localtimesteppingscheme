#include <stdio.h>
#include <stdlib.h>
#include "computeFlux.h"
#include "macroDefinition.h"
#include "updateFlux.h"
#include "visual.h"

int main(int argc, char **argv)
{
	/* initial value set up */
	int cellsize = 1;
	int n_grid = 50;
	int length = n_grid * cellsize;
	int totalNumberofTimeStep = 1000;
	int plottingStep = 1;
	double dt = 0.1;
	double dt_dx = dt/cellsize;

	const char *szProblem;
	szProblem = "result";

	/* Initialisation & memory allocation */
	double amax;
	double *h, *u, *v, *F, *G, *U;

	h = malloc(n_grid*n_grid*sizeof(double));
	u = malloc(n_grid*n_grid*sizeof(double));
	v = malloc(n_grid*n_grid*sizeof(double));
	F = malloc((n_grid+1)*n_grid*3*sizeof(double));
	G = malloc((n_grid+1)*n_grid*3*sizeof(double));
	U = malloc(n_grid*n_grid*3*sizeof(double));

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
		}
	}
	/*initialise h*/
	for (int x = 0; x < 10; ++x)
	{
		for (int y = n_grid - 20; y < n_grid-10; ++y)
		{
			h[x*n_grid + y] = 1.0;
			U[ (x*n_grid + y)*3] = h[x*n_grid + y];
		}
	}

	for (int x = 0; x < n_grid + 1; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
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
			for (int i = 0; i < 3; ++i)
			{
				G[ (x*(n_grid+1) + y)*3 + i ] = 0.0;
			}
		}
	}

	for (int i = 0; i < totalNumberofTimeStep; ++i)
	{
	// for (int i = 0; i < 3; ++i)
	// {
	// 	for (int x = 0; x < n_grid; ++x)
	// 	{
	// 		for (int y = 0; y < n_grid; ++y)
	// 		{
				
	// 				printf("U: %lf \t", U[(x*n_grid + y)*3 + i]);
				
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }
		if(i % plottingStep == 0){
        	write_vtkFile(szProblem, i, length, n_grid, n_grid, cellsize, cellsize, U);
        }

		// tile creation

		/* compute fluxes*/
		computeFlux(U, F, G, n_grid, &amax);

		/* updating the fluxes*/
		updateFlux(U, F, G, n_grid, dt_dx);

		//printf("Time Step = %d, amax = %lf \n", i, amax);
		printf("Time Step = %d, Courant Number = %lf \n", i, amax * dt_dx* 2 );
		/* write vtk file*/

	}

	/* memory deallocation */
	free(h);
	free(u);
	free(v);
	free(F);
	free(G);
	free(U);
}
