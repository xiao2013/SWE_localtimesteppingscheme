#include "computeTimeStep.h"

void calculateTimeStep(double *lambdaf, double *lambdag, double *dt, int *levelc,
	int *levelf, int *levelg, int LTS_levels, int cellsize, int n_grid, double dt0, double *maxCrnumber)
{
	double crmax = 0.;
	double crdum, crtest;

	double *cr;

	cr = malloc(n_grid*n_grid*sizeof(double));

	for (int x = 0; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			cr[x*n_grid + y] = 0.0;
		}
	}

	if (LTS_levels > 1)
	{
		/* Compute Courant number and assign level id to cells */
		for (int x = 0; x < n_grid; ++x)
		{
			for (int y = 0; y < n_grid; ++y)
			{
				crdum = dt0*(fmax(lambdaf[x*n_grid + y],lambdaf[(x+1)*n_grid + y])
				+ fmax(lambdag[x*(n_grid+1) + y],lambdag[x*(n_grid+1) + y + 1]))/cellsize;
				// printf("crdum[%d][%d]: %lf\t", x, y, crdum);
				for (int i = 1; i <= LTS_levels - 1; ++i)
				{
					crtest=0.8*pow(2,-i);
					// printf("crtest: %lf\n", crtest);
					if (crdum >= crtest)
					{
						// printf("ok!1\n");
						levelc[x*n_grid + y] = i;
						dt[x*n_grid + y] = dt0*pow(2, i - 1);
						cr[x*n_grid + y] = crdum*dt[x*n_grid + y]/dt0; // @dt0: dt is "adapted" here
											 // printf("levelc[%d][%d]: %d\t", x, y, levelc[x*n_grid + y]);

						break;
					}
					levelc[x*n_grid + y] = LTS_levels;
					dt[x*n_grid + y] = dt0*pow(2, LTS_levels - 1); //@dt0: dt is "adapted" here
					// printf("dt[%d][%d]: %lf\t", x, y, dt[x*n_grid + y]);
					cr[x*n_grid + y] = crdum*dt[x*n_grid + y]/dt0;

				}
				// printf("\n");

				if (cr[x*n_grid + y]>crmax)
				{
					crmax=cr[x*n_grid + y];
				}
			}
		}

		/* Assign a level id to edges of mesh */
		for (int x = 1; x < n_grid; ++x)
		{
			for (int y = 0; y < n_grid; ++y)
			{
				levelf[x*n_grid + y] = fmin(levelc[(x-1)*n_grid + y], levelc[x*n_grid + y]);
			}
		}
		for (int y = 0; y < n_grid; ++y)
		{
			levelf[y]=levelc[y];
			levelf[n_grid*n_grid + y] = levelc[(n_grid-1)*n_grid + y];
		}

		for (int y = 1; y < n_grid; ++y)
		{
			for (int x = 0; x < n_grid; ++x)
			{
				levelg[x*(n_grid+1) + y] = fmin(levelc[x*n_grid + y-1], levelc[x*n_grid + y]);
			}
		}
		for (int x = 0; x < n_grid; ++x)
		{
			levelg[x*(n_grid+1)] = levelc[x*n_grid];
			levelg[x*(n_grid+1) + n_grid] = levelc[x*n_grid + n_grid-1];
		}

		/* Create a layer of buffer cells */
		for (int x = 0; x < n_grid; ++x)
		{
			for (int y = 0; y < n_grid; ++y)
			{
				levelc[x*n_grid + y] = min4(levelf[x*n_grid + y] , levelf[(x+1)*n_grid + y],
					levelg[x*(n_grid+1) + y],levelg[x*(n_grid+1) + y+1]);
				dt[x*n_grid + y] = dt0*pow(2, levelc[x*n_grid + y]-1); // @dt0: dt is "adapted" here
				// printf("dt[%d][%d]: %lf\t", x, y, dt[x*n_grid + y]);

			}
			// printf("\n");
		}
	}
	else
	{
		for (int x = 0; x < n_grid; ++x)
		{
			for (int y = 0; y < n_grid; ++y)
			{
				cr[x*n_grid + y] = dt0*(fmax(lambdaf[x*n_grid + y],lambdaf[(x+1)*n_grid + y])
				+ fmax(lambdag[x*(n_grid+1) + y],lambdag[x*(n_grid+1) + y + 1]))/cellsize;
				if (cr[x*n_grid + y] > crmax)
				{
					crmax=cr[x*n_grid + y];
				}
			}
		}
	}
	*maxCrnumber = crmax;
	free(cr);
}

double min4(double a, double b, double c, double d)
{
	return fmin(a, fmin(b, fmin(c, d)));
}