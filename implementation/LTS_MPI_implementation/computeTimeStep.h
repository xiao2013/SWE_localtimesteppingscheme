#ifndef __COMPUTETIMESTEP_H__
#define __COMPUTETIMESTEP_H__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include "macroDefinition.h"
#include "helper.h"

void calculateTimeStep(double *lambdaf, double *lambdag, double *dt, int *levelc,
	int *levelf, int *levelg, int LTS_levels, int cellsize, int n_grid, double dt0, double *crmax);

double min4(double a, double b, double c, double d);

#endif