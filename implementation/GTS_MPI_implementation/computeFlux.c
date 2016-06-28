#include "computeFlux.h"
#include "computeRoeSolver.h"
#include <math.h>
#include <stdio.h>

int currentIndex(int x, int y, int xlength){
	return x * xlength + y;
}

void computeFlux(double *U, double *F, double *G, int n_grid, double *maxValue, int mex, int mey, int npx)
{
	/* compute fluxes in x direction*/
	int g_n = 0;
	int f_n = 1;
	double amax = 0;
	double a;
	double hl, hr, ul, ur, vl, vr;
	double local_F[3]= {0.};
	int x, y;
	if (mex != 0)
	{
		for (int x = 1; x < n_grid+1; ++x)
		{
			for (int y = 1; y < n_grid+1; ++y)
			{
				hl = U[currentIndex(x-1, y, (n_grid+2)) * 3];
				hr = U[currentIndex(x, y, (n_grid+2)) * 3];
				ul = U[currentIndex(x-1, y, (n_grid+2)) * 3 + 1] / hl;
				ur = U[currentIndex(x, y, (n_grid+2)) * 3 + 1] / hr;
				vl = U[currentIndex(x-1, y, (n_grid+2)) * 3 + 2] / hl;
				vr = U[currentIndex(x, y, (n_grid+2)) * 3 + 2] / hr;
				// printf("hl[%d][%d]: %lf\n",x, y, hl);
				// printf("hr[%d][%d]: %lf\n",x, y, hr);
				// printf("ul[%d][%d]: %lf\n",x, y, ul);
				// printf("ur[%d][%d]: %lf\n",x, y, ur);
				// printf("vl[%d][%d]: %lf\n",x, y, vl);
				// printf("vr[%d][%d]: %lf\n",x, y, vr);

				/* compute roe solver */
				//printf("x: %d, y : %d\n",x,y );
				computeRoeSolver(hl,hr,ul,ur,vl,vr,g_n,f_n, local_F, &a);
				//printf("a: %lf\n", a);

				for (int i = 0; i < 3; ++i)
				{
					F[(x*(n_grid+2) + y)*3 + i] = local_F[i];
				}
				amax = fmax(a, amax);

			}
		}
	}
	else // mex == 0
	{
		for (int x = 2; x <= n_grid+1; ++x)
		{
			for (int y = 1; y < n_grid+1; ++y)
			{
				hl = U[currentIndex(x-1, y, (n_grid+2)) * 3];
				hr = U[currentIndex(x, y, (n_grid+2)) * 3];
				ul = U[currentIndex(x-1, y, (n_grid+2)) * 3 + 1] / hl;
				ur = U[currentIndex(x, y, (n_grid+2)) * 3 + 1] / hr;
				vl = U[currentIndex(x-1, y, (n_grid+2)) * 3 + 2] / hl;
				vr = U[currentIndex(x, y, (n_grid+2)) * 3 + 2] / hr;
				// printf("hl[%d][%d]: %lf\n",x, y, hl);
				// printf("hr[%d][%d]: %lf\n",x, y, hr);
				// printf("ul[%d][%d]: %lf\n",x, y, ul);
				// printf("ur[%d][%d]: %lf\n",x, y, ur);
				// printf("vl[%d][%d]: %lf\n",x, y, vl);
				// printf("vr[%d][%d]: %lf\n",x, y, vr);

				/* compute roe solver */
				//printf("x: %d, y : %d\n",x,y );
				computeRoeSolver(hl,hr,ul,ur,vl,vr,g_n,f_n, local_F, &a);
				//printf("a: %lf\n", a);

				for (int i = 0; i < 3; ++i)
				{
					F[(x*(n_grid+2) + y)*3 + i] = local_F[i];
				}
				amax = fmax(a, amax);

			}
		}
	}
	
	// for (int i = 0; i < 3; ++i)
	// {
	// 	for (int x = 0; x < n_grid + 1; ++x)
	// 	{
	// 		for (int y = 0; y < n_grid; ++y)
	// 		{
	// 			printf("F:[%d][%d]: %lf\t", x, y,F[(x*n_grid + y)*3 + i]);
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }
	// printf("part 1 done! amax: %lf\n", amax);
	if (mex == 0)
	{
		/* Wall boundary on x = 0 */
		x = 1;
		for (int y = 1; y < n_grid+1; ++y)
		{
			hr = U[currentIndex(x, y, n_grid+2) * 3];
			ur = U[currentIndex(x, y, n_grid+2) * 3 + 1] / hr;
			vr = U[currentIndex(x, y, n_grid+2) * 3 + 2] / hr;
			/* compute roe solver */
			computeRoeSolver(hr,hr,-ur,ur,vr,vr,g_n,f_n, local_F, &a);
			for (int i = 0; i < 3; ++i)
			{
				F[(x*(n_grid+2) + y)*3 + i] = local_F[i];
			}
			amax = fmax(a, amax);
		}
	}

	//printf("x = 0 done! amax: %lf\n", amax);


	/* Wall boundary on x = n_grid -1 */
	if (mex == (npx - 1))
	{
		for (int y = 1; y < n_grid+1; ++y)
		{
			hl = U[currentIndex(n_grid, y, n_grid+2) * 3];
			ul = U[currentIndex(n_grid, y, n_grid+2) * 3 + 1] / hl;
			vl = U[currentIndex(n_grid, y, n_grid+2) * 3 + 2] / hl;
			/* compute roe solver */
			computeRoeSolver(hl,hl,ul,-ul,vl,vl,g_n,f_n, local_F, &a);
			for (int i = 0; i < 3; ++i)
			{
				F[((n_grid+1)*(n_grid+2) + y)*3 + i] = local_F[i];
			}
			amax = fmax(a, amax);
		}
	}
	// printf("x = n_grid done! amax: %lf\n", amax);

		
	/* compute fluxes in y direction */
	g_n = 1;
	f_n = 0;
	if (mey == 0)
	{

		for (int y = 2; y <= n_grid+1; ++y)
		{
			// @reminder: MPI the element indexed by the current y value (y=0) was already present from initiation of the tiles in the buffer column
			for (int x = 1; x < n_grid+1; ++x)
			{

				hl = U[currentIndex(x, y - 1, n_grid+2) * 3];
				hr = U[currentIndex(x, y, n_grid+2) * 3];
				//printf("hl[%d][%d]: %lf\n",x, y, hl);
				ul = U[currentIndex(x, y - 1, n_grid+2) * 3 + 1] / hl;
				ur = U[currentIndex(x, y, n_grid+2) * 3 + 1] / hr;
				vl = U[currentIndex(x, y - 1, n_grid+2) * 3 + 2] / hl;
				vr = U[currentIndex(x, y, n_grid+2) * 3 + 2] / hr;
				/* compute roe solver */
				computeRoeSolver(hl,hr,ul,ur,vl,vr,g_n,f_n, local_F, &a);
				for (int i = 0; i < 3; ++i)
				{
					G[(x*(n_grid+2) + y)*3 + i] = local_F[i];
				}
				amax = fmax(a, amax);
			}

			// @reminder: copy
		}
	}
	else // mey != 0
	{
		for (int y = 1; y < n_grid+1; ++y)
		{
			// @reminder: MPI the element indexed by the current y value (y=0) was already present from initiation of the tiles in the buffer column
			for (int x = 1; x < n_grid+1; ++x)
			{

				hl = U[currentIndex(x, y - 1, n_grid+2) * 3];
				hr = U[currentIndex(x, y, n_grid+2) * 3];
				//printf("hl[%d][%d]: %lf\n",x, y, hl);
				ul = U[currentIndex(x, y - 1, n_grid+2) * 3 + 1] / hl;
				ur = U[currentIndex(x, y, n_grid+2) * 3 + 1] / hr;
				vl = U[currentIndex(x, y - 1, n_grid+2) * 3 + 2] / hl;
				vr = U[currentIndex(x, y, n_grid+2) * 3 + 2] / hr;
				/* compute roe solver */
				computeRoeSolver(hl,hr,ul,ur,vl,vr,g_n,f_n, local_F, &a);
				for (int i = 0; i < 3; ++i)
				{
					G[(x*(n_grid+2) + y)*3 + i] = local_F[i];
				}
				amax = fmax(a, amax);
			}

			// @reminder: copy
		}
	}
	// printf("part 2 done! amax: %lf\n", amax);
	if (mey == 0)
	{
		/* Wall boundary on y = 0 */
		y = 1;
		for (int x = 1; x < n_grid+1; ++x)
		{
			hr = U[currentIndex(x, y, n_grid+2) * 3];
			ur = U[currentIndex(x, y, n_grid+2) * 3 + 1] / hr;
			vr = U[currentIndex(x, y, n_grid+2) * 3 + 2] / hr;
			/* compute roe solver */
			computeRoeSolver(hr,hr,ur,ur,-vr,vr,g_n,f_n, local_F, &a);
			for (int i = 0; i < 3; ++i)
			{
				G[(x*(n_grid+2) + y)*3 + i] = local_F[i];
			}
			amax = fmax(a, amax);
		}
	}

	//printf("\n");

	// printf("y = 0 done! amax: %lf\n", amax);
	if (mey == (npx -1))
	{
		y = n_grid;
		/* Wall boundary on y = n_grid -1*/
		for (int x = 1; x < n_grid+1; ++x)
		{
			hl = U[currentIndex(x, y, n_grid+2) * 3];
			ul = U[currentIndex(x, y, n_grid+2) * 3 + 1] / hl;
			vl = U[currentIndex(x, y, n_grid+2) * 3 + 2] / hl;
			/* compute roe solver */
			computeRoeSolver(hl,hl,ul,ul,vl,-vl,g_n,f_n, local_F, &a);
			for (int i = 0; i < 3; ++i)
			{
				G[(x*(n_grid+2) + y+1)*3 + i] = local_F[i];
			}
			amax = fmax(a, amax);
		}
	}

	// for (int i = 0; i < 3; ++i)
	// {
	// 	for (int x = 0; x < n_grid; ++x)
	// 	{
	// 		for (int y = 0; y < n_grid + 1; ++y)
	// 		{
	// 			printf("G:[%d][%d]: %lf\t", x, y,G[(x*(n_grid + 1) + y)*3 + i]);			
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }
	// printf("y = n_grid -1 done! amax: %lf\n", amax);
	*maxValue = amax;
	//printf("amax :%lf\n", amax );
}
