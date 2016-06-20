#include "computeFlux.h"
#include "computeRoeSolver.h"
#include <math.h>
#include <stdio.h>

void printMatrix(int n_grid, int *matrix)
{
	for (int x = 0; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			printf("matrix[%d][%d]: %d\t",x,y, matrix[x*n_grid + y] );
		}
		printf("\n");
	}
}


int currentIndex(int x, int y, int xlength){
	return x * xlength + y;
}

void computeFlux(double *U, double *F, double *G, int *levelf, int *levelg, double *lambdaf,
				 double *lambdag,int n_grid, int level)
{
	/* compute fluxes in x direction*/
	int g_n = 0;
	int f_n = 1;
	// double amax = 0;
	// double a;
	double hl, hr, ul, ur, vl, vr;
	double local_F[3]= {0.};

	// printMatrix(n_grid, levelf);

	for (int x = 1; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			if (levelf[x*n_grid + y] <= level)
			{
				hl = U[currentIndex(x-1, y, n_grid) * 3];
				hr = U[currentIndex(x, y, n_grid) * 3];
				ul = U[currentIndex(x-1, y, n_grid) * 3 + 1] / hl;
				ur = U[currentIndex(x, y, n_grid) * 3 + 1] / hr;
				vl = U[currentIndex(x-1, y, n_grid) * 3 + 2] / hl;
				vr = U[currentIndex(x, y, n_grid) * 3 + 2] / hr;
				//if (x == 1 && y==4)
				// {
				// 	printf("hl[%d][%d]: %lf\n",x, y, hl);
				// printf("hr[%d][%d]: %lf\n",x, y, hr);
				// printf("ul[%d][%d]: %lf\n",x, y, ul);
				// printf("ur[%d][%d]: %lf\n",x, y, ur);
				// printf("vl[%d][%d]: %lf\n",x, y, vl);
				// printf("vr[%d][%d]: %lf\n",x, y, vr);
				
				/* compute roe solver */
				//printf("x: %d, y : %d\n",x,y );
				computeRoeSolver(hl,hr,ul,ur,vl,vr,g_n,f_n, local_F, &lambdaf[x*n_grid+y]);
				
				for (int i = 0; i < 3; ++i)
				{
					F[(x*n_grid + y)*3 + i] = local_F[i];
					// printf("local_F: %lf\n", local_F[i]);
				}
				// }
				//amax = fmax(a, amax);
			}
		}
	}
	// printf("F[1][4][1]:%lf\n", F[(1*n_grid + 4)*3 + 1]);

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

	/* Wall boundary on x = 0 */
	for (int y = 0; y < n_grid; ++y)
	{
		if (levelf[0*n_grid + y] <= level)
		{
			hr = U[currentIndex(0, y, n_grid) * 3];
			ur = U[currentIndex(0, y, n_grid) * 3 + 1] / hr;
			vr = U[currentIndex(0, y, n_grid) * 3 + 2] / hr;
			/* compute roe solver */
			computeRoeSolver(hr,hr,-ur,ur,vr,vr,g_n,f_n, local_F, &lambdaf[0*n_grid+y]);
			for (int i = 0; i < 3; ++i)
			{
				F[(0*n_grid + y)*3 + i] = local_F[i];
			}
			//amax = fmax(a, amax);
		}
	}
	//printf("x = 0 done! amax: %lf\n", amax);


	/* Wall boundary on x = n_grid -1 */
	for (int y = 0; y < n_grid; ++y)
	{
		if (levelf[n_grid*n_grid + y] <= level)
		{
			hl = U[currentIndex(n_grid - 1, y, n_grid) * 3];
			ul = U[currentIndex(n_grid - 1, y, n_grid) * 3 + 1] / hl;
			vl = U[currentIndex(n_grid - 1, y, n_grid) * 3 + 2] / hl;
			/* compute roe solver */
			computeRoeSolver(hl,hl,ul,-ul,vl,vl,g_n,f_n, local_F, &lambdaf[n_grid*n_grid+y]);
			for (int i = 0; i < 3; ++i)
			{
				F[(n_grid*n_grid + y)*3 + i] = local_F[i];
			}
			//amax = fmax(a, amax);
		}
	}
	// printf("x = n_grid done! amax: %lf\n", amax);

		
	/* compute fluxes in y direction */
	g_n = 1;
	f_n = 0;
	for (int y = 1; y < n_grid; ++y)
	{
		for (int x = 0; x < n_grid; ++x)
		{
			if (levelg[x*(n_grid+1) + y] <= level)
			{
				// printf("ok!\n");
				hl = U[currentIndex(x, y - 1, n_grid) * 3];
				hr = U[currentIndex(x, y, n_grid) * 3];
				//printf("hl[%d][%d]: %lf\n",x, y, hl);
				ul = U[currentIndex(x, y - 1, n_grid) * 3 + 1] / hl;
				// printf("ul[%d][%d]:%lf\n",x,y,ul );
				ur = U[currentIndex(x, y, n_grid) * 3 + 1] / hr;
				vl = U[currentIndex(x, y - 1, n_grid) * 3 + 2] / hl;
				vr = U[currentIndex(x, y, n_grid) * 3 + 2] / hr;
				/* compute roe solver */
				computeRoeSolver(hl,hr,ul,ur,vl,vr,g_n,f_n, local_F, &lambdag[x*(n_grid+1)+y]);
				for (int i = 0; i < 3; ++i)
				{
					G[(x*(n_grid+1) + y)*3 + i] = local_F[i];
				}
				//amax = fmax(a, amax);
			}
		}
	}

	// printf("part 2 done! amax: %lf\n", amax);

	/* Wall boundary on y = 0 */
	for (int x = 0; x < n_grid; ++x)
	{
		if (levelg[x*(n_grid+1) + 0] <= level)
		{
			hr = U[currentIndex(x, 0, n_grid) * 3];
			ur = U[currentIndex(x, 0, n_grid) * 3 + 1] / hr;
			vr = U[currentIndex(x, 0, n_grid) * 3 + 2] / hr;
			/* compute roe solver */
			computeRoeSolver(hr,hr,ur,ur,-vr,vr,g_n,f_n, local_F, &lambdag[x*(n_grid+1)+0]);
			for (int i = 0; i < 3; ++i)
			{
				G[(x*(n_grid+1) + 0)*3 + i] = local_F[i];
			}
			//amax = fmax(a, amax);
		}
	}
	// printf("\n");

	// printf("y = 0 done! amax: %lf\n", amax);

	/* Wall boundary on y = n_grid -1*/
	for (int x = 0; x < n_grid; ++x)
	{
		if (levelg[x*(n_grid+1) + n_grid] <= level)
		{
			hl = U[currentIndex(x, n_grid - 1, n_grid) * 3];
			ul = U[currentIndex(x, n_grid - 1, n_grid) * 3 + 1] / hl;
			vl = U[currentIndex(x, n_grid - 1, n_grid) * 3 + 2] / hl;
			/* compute roe solver */
			computeRoeSolver(hl,hl,ul,ul,vl,-vl,g_n,f_n, local_F, &lambdag[x*(n_grid+1)+n_grid]);
			for (int i = 0; i < 3; ++i)
			{
				G[(x*(n_grid+1) + n_grid)*3 + i] = local_F[i];
			}
			//amax = fmax(a, amax);
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
	// *maxValue = amax;
	//printf("amax :%lf\n", amax );
}
