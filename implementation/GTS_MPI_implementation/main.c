#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "computeFlux.h"
#include "macroDefinition.h"
#include "updateFlux.h"
#include "visual.h"

# define physx( i, mex, NX ) ( (i) + mex*NX )
# define physy( j, mey, NY ) ( (j) + mey*NY )

int main(int argc, char **argv)
{

    /* initialise MPI */
    int size, rank;
    int mex, mey;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* initial value set up */
    int cellsize = 1;
    int n_grid = 40;
    int length = n_grid * cellsize;
    int totalNumberofTimeStep = 1;
    int plottingStep = 1;
    double dt = 0.1;
    double dt_dx = dt/cellsize;
    const char *szProblem;
    szProblem = "result";

    /* for MPI tiles */
    int npx = sqrt(size);    // for making parallel processed tiles, npx is the square length of the tile grid e.g. if there are 9 processors then npx = 3)
    int l_grid = n_grid/npx; // l_grid = local grid, e.g. if n_grid = 50, and npx = 5, then the number of global grid elements per tile is 10 (in one dimension)

    // npx = 4/4 = 1

    printf("npx: %i\tn_grid: %i\tl_grid: %i\n", npx, n_grid,l_grid); // @ANDY:debug
    /* for mpi_subarray */
    int fullsize[2] = {n_grid, n_grid};
    int localsize[2] = {l_grid, l_grid};
    int starts[2] = {0, 0};
    /* create 2D mpi processor virtual topology */
    MPI_Comm mpi_comm_cart, mpi_comm_x, mpi_comm_y;
    int dims[2];
    int periods[2];
    int coords[2];
    int reorder = 0;
    dims[0] = npx; /* number of processor in x direction */
    dims[1] = npx; /* number of processor in y direction */
    periods[0] = 0;
    periods[1] = 0;

    MPI_Cart_create (MPI_COMM_WORLD, 2, dims, periods, reorder, &mpi_comm_cart);
    /* Build the sub-communicators along X and Y */
    coords[0] = 1;
    coords[1] = 0;
    MPI_Cart_sub (mpi_comm_cart, coords, &mpi_comm_x);
    coords[0] = 0;
    coords[1] = 1;
    MPI_Cart_sub (mpi_comm_cart, coords, &mpi_comm_y);
    /* Rank along X, Y directions */
    MPI_Comm_rank (mpi_comm_x, &mex);
    MPI_Comm_rank (mpi_comm_y, &mey);

    /* Initialisation & memory allocation */
    double amax;
    double *h, *u, *v, *F, *G, *U, *U_global, *sendArray;
    int xmax, xmin, ymin, ymax;

    xmin = 0;
    xmax = 10;
    ymin = n_grid - 20;
    ymax = n_grid - 10;
    
    // h = malloc(n_grid*n_grid*sizeof(double));
    // u = malloc(n_grid*n_grid*sizeof(double));
    // v = malloc(n_grid*n_grid*sizeof(double));
    // F = malloc((n_grid+1)*n_grid*3*sizeof(double));
    // G = malloc((n_grid+1)*n_grid*3*sizeof(double));
    U_global = malloc(n_grid*n_grid*3*sizeof(double));
    h = malloc((l_grid+2)*(l_grid+2)*sizeof(double));
    u = malloc((l_grid+2)*(l_grid+2)*sizeof(double));
    v = malloc((l_grid+2)*(l_grid+2)*sizeof(double));
    F = malloc((l_grid+2)*(l_grid+2)*3*sizeof(double));
    G = malloc((l_grid+2)*(l_grid+2)*3*sizeof(double));
    U = malloc((l_grid+2)*(l_grid+2)*3*sizeof(double));
    sendArray = malloc(l_grid*l_grid*3*sizeof(double));

<<<<<<< HEAD
    //printf("mex: %d, mey: %d\n rank: %d\n", mex, mey, rank);  // @ANDY:debug
    //printf("l_grid: %d\n", l_grid);                           // @ANDY:debug
=======
    printf("mex: %d, mey: %d rank: %d\n", mex, mey, rank);
    printf("x: %d, y: %d\n", physx(0, mex, l_grid), physy(0, mey, l_grid));
    printf("l_grid: %d\n", l_grid);
>>>>>>> f3e4983b442065f4e165c5844a69f6e6cc6379da
    for (int x = 1; x < l_grid+1; ++x) // so x = 0 and x = l_grid+2 are unallocated "ghost layers" 
    {
        for (int y = 1; y < l_grid+1; ++y)
        {
            h[x*(l_grid+2) + y] = 0.1;
            u[x*(l_grid+2) + y] = 0.0;
            v[x*(l_grid+2) + y] = 0.0;
            U[ (x*(l_grid+2) + y)*3]     = h[x*(l_grid+2) + y];
            U[ (x*(l_grid+2) + y)*3 + 1] = u[x*(l_grid+2) + y] * h[x*(l_grid+2) + y];
            U[ (x*(l_grid+2) + y)*3 + 2] = v[x*(l_grid+2) + y] * h[x*(l_grid+2) + y];
            /* initialise the shock wave by lifting up parts of the water */
            if (physx(x, mex, (l_grid)) < (xmax+1) 
                && physx(x, mex, (l_grid)) > (xmin) 
                && physy(y, mey, (l_grid)) < (ymax+1)
                && physy(y, mey, (l_grid)) > (ymin))
            // if (x < (xmax+1) 
            //     && x > (xmin) 
            //     && y < 41
            //     && y > 30)
            {
                // counter++;
                h[x*(l_grid + 2) + y] = 1.0;
                U[ (x*(l_grid + 2) + y)*3] = h[x*(l_grid + 2) + y];
            }
        }
    }
    // printf("counter: %d\n", counter);
    /*initialise h*/

    // for (int x = 1; x < 11; ++x)
    // {
    //     for (int y = n_grid - 19; y < n_grid-9; ++y)
    //     {
    //         h[x*(l_grid+2) + y] = 1.0;
    //         U[ (x*(l_grid+2) + y)*3] = h[x*(l_grid+2) + y];
    //     }
    // }

    for (int x = 1; x < l_grid+2; ++x)
    {
        for (int y = 1; y < l_grid+1; ++y)
        {
            for (int i = 0; i < 3; ++i)
            {
                F[ (x*(l_grid + 2) + y)*3 + i ] = 0.0;
            }
        }
    }
    for (int x = 1; x < l_grid+1; ++x)
    {
        for (int y = 1; y < l_grid+2; ++y)
        {
            for (int i = 0; i < 3; ++i)
            {
                G[ (x*(l_grid + 2) + y)*3 + i ] = 0.0;
            }
        }
    }

    for (int i = 0; i < totalNumberofTimeStep; ++i)
    {
    // for (int i = 0; i < 3; ++i)
    // {
    //  for (int x = 0; x < n_grid; ++x)
    //  {
    //      for (int y = 0; y < n_grid; ++y)
    //      {
                
    //              printf("U: %lf \t", U[(x*n_grid + y)*3 + i]);
                
    //      }
    //      printf("\n");
    //  }
    //  printf("\n");
    // }
        // if(i == 999){
        //         for (int x = 0; x < n_grid; ++x)
        //         {
        //             for (int y = 0; y < n_grid; ++y)
        //             {
        //                 printf("%lf\t",U[ (x*n_grid + y)*3]);
        //             }
        //             printf("\n");
        //         }
        // }

        // //tile creation

        // /* compute fluxes*/
        // computeFlux(U, F, G, n_grid, &amax);



        // /* updating the fluxes*/
        // updateFlux(U, F, G, n_grid, dt_dx);

        //printf("Time Step = %d, amax = %lf \n", i, amax);
        //printf("Time Step = %d, Courant Number = %lf \n", i, amax * dt_dx* 2 );
        /* gather all information from all processors to main */
        int counter = 0;
       
        if (i % plottingStep == 0)
        {
            for (int x = 0; x < l_grid; ++x) // so x = 0 and x = l_grid+2 are unallocated "ghost layers" 
            {
                for (int y = 0; y < l_grid; ++y)
                {
                    if (U[ ((x+1)*(l_grid + 2) + (y+1))*3 ] == 1)
                    {
                        // printf("%d\n", rank);s
                        counter++;
                    }
                    sendArray[ (x*l_grid + y)*3 ]     = U[ ((x+1)*(l_grid + 2) + (y+1))*3 ];
                    sendArray[ (x*l_grid + y)*3 + 1 ] = U[ ((x+1)*(l_grid + 2) + (y+1))*3 + 1 ];
                    sendArray[ (x*l_grid + y)*3 + 2 ] = U[ ((x+1)*(l_grid + 2) + (y+1))*3 + 2 ];
                }
            }
            printf("rank: %d counter:%d\n", rank, counter);
            MPI_Gather( sendArray, l_grid*l_grid*sizeof(double)*3, MPI_DOUBLE, U_global, l_grid*l_grid*sizeof(double)*3, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
            /* write vtk file*/
            if (rank == 0)
            {
                int counter_rank0 = 0;
                for (int x = 0; x < n_grid; ++x) // so x = 0 and x = l_grid+2 are unallocated "ghost layers" 
                {
                    for (int y = 0; y < n_grid; ++y)
                    {
                        for (int i = 0; i < 3; ++i)
                        {
                            if (U_global[ ((x)*(n_grid) + (y))*3 ] == 1)
                            {
                                counter_rank0++;
                            }
                        }


                    }
                }
                printf("%d\n", counter_rank0);
            }

            // if( rank == 0){
            //     printf("ok!\n");
            //     write_vtkFile(szProblem, i, length, n_grid, n_grid, cellsize, cellsize, U_global);
            // }
        }
        printf("rank number: %i\n",rank);
    }

    /* memory deallocation */
    free(h);
    free(u);
    free(v);
    free(F);
    free(G);
    free(U);
    free(U_global);
    free(sendArray);

    MPI_Finalize();
}
