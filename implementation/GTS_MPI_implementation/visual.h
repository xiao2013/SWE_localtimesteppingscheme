#ifndef __VISUAL_H__
#define __VISUAL_H__

void write_vtkFile( const char *szProblem,int timeStepNumber, int gridsize, int imax, int jmax,
                    double dx, double dy, double *U);

void write_vtkHeader( FILE *fp, int imax, int jmax, double dx, double dy);

void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, double dx, double dy);

#endif