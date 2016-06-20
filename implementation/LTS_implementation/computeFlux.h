#ifndef __COMPUTEFLUX_H__
#define __COMPUTEFLUX_H__

void computeFlux(double *U, double *F, double *G, int *levelf, int *levelg, double *lambdaf,
	double *lambdag, int n_grid, int level);

#endif