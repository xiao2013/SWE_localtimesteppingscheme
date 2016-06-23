#ifndef __COMPUTEROESOLVER_H__
#define __COMPUTEROESOLVER_H__

void computeRoeSolver(double hl, double hr, double ul, double ur, double vl, double vr,
					 int g_n, int f_n, double *F, double *amax);

void computeRes(double R[3][3], double A[3][3], double *dW, double *res);

#endif