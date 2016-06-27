void updateFlux(double *U, double *F, double *G, int n_grid, double dt_dx)
{
	for (int x = 1; x < n_grid+1; ++x)
	{
		for (int y = 1; y < n_grid+1; ++y)
		{
			for (int i = 0; i < 3; ++i)
			{
				U[ (x*(n_grid+2) + y)*3 + i] = U[ (x*(n_grid+2) + y)*3 + i] - dt_dx* (F[ ((x+1)*(n_grid+2) + y)*3 + i ] 
					- F[ (x*(n_grid+2) + y)*3 + i ] + G[ (x*(n_grid+2) + y + 1)*3 + i ] - G[ (x*(n_grid+2) + y)*3 + i ]);
				//                                            ^(n_grid+1)
			}
		}
	}
}