#include <stdio.h>

void updateFlux(double *U, double *F, double *G, int *levelc, int n_grid, double *dt, int level)
{
		// 	for (int x = 0; x < n_grid; ++x)
		// {
		// for (int y = 0; y < n_grid; ++y)
		// {
		// 	printf("levelc[%d][%d]: %d\t",x,y, levelc[x*n_grid + y] );
		// }
		// printf("\n");
		// }
		// 	for (int x = 0; x < n_grid; ++x)
		// {
		// for (int y = 0; y < n_grid; ++y)
		// {
		// 	printf("dt[%d][%d]: %lf\t",x,y, dt[x*n_grid + y] );
		// }
		// printf("\n");
		// }
	// 	for (int i = 0; i < 3; ++i)
	// {
	// 	for (int x = 0; x < n_grid + 1; ++x)
	// 	{
	// 		for (int y = 0; y < n_grid; ++y)
	// 		{
	// 			printf("F[%d][%d]: %lf \t",x,y, F[(x*n_grid + y)*3 + i]);
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }

	// 		for (int i = 0; i < 3; ++i)
	// {
	// 	for (int x = 0; x < n_grid; ++x)
	// 	{
	// 		for (int y = 0; y < n_grid + 1; ++y)
	// 		{
	// 			printf("G: %lf \t", G[(x*(n_grid+1) + y)*3 + i]);
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }

	for (int x = 0; x < n_grid; ++x)
	{
		for (int y = 0; y < n_grid; ++y)
		{
			if (levelc[x*n_grid + y] <= level) // @ANDY: if the current grid point needs to be updated (i.e. sufficiently small local dt) then it will be updated. Grid points that have large dt values would not need to be calculated as frequently as points with small local dt. 
			{
				// printf("[%d][%d]\n", x, y);
				for (int i = 0; i < 3; ++i)
				{
					U[ (x*n_grid + y)*3 + i] = U[ (x*n_grid + y)*3 + i] - dt[x*n_grid + y]* (F[ ((x+1)*n_grid + y)*3 + i ] 
						- F[ (x*n_grid + y)*3 + i ] + G[ (x*(n_grid+1) + y + 1)*3 + i ] - G[ (x*(n_grid+1) + y)*3 + i ]);
				}
			}
		}
	}
}