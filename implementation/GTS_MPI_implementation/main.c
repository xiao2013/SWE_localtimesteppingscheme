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

    /* initialise */
    int debug = 0;  //@ANDY:debug @debug
    int x, y;
    int size, rank;
    int mex, mey;  // (0,0) in (x,y), is at the bottom left of the grid of tiles, whre tiles are processors
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);   

    /* initial value set up */
    int cellsize = 1;
    int n_grid = 4;
    int length = n_grid * cellsize;
    int totalNumberofTimeStep = 1000;
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
    int localsize[2] = {l_grid , l_grid };
    int starts[2] = {0, 0};
    MPI_Datatype subarrayType, type;

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
    /* type for global array displacement */
    MPI_Type_create_subarray(2, fullsize, localsize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    MPI_Type_create_resized(type, 0, l_grid*sizeof(double), &subarrayType);
    MPI_Type_commit(&subarrayType);
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
    U_global = malloc(n_grid*n_grid*sizeof(double));
    h = malloc((l_grid+2)*(l_grid+2)*sizeof(double));
    u = malloc((l_grid+2)*(l_grid+2)*sizeof(double));
    v = malloc((l_grid+2)*(l_grid+2)*sizeof(double));
    F = malloc((l_grid+2)*(l_grid+2)*3*sizeof(double));
    G = malloc((l_grid+2)*(l_grid+2)*3*sizeof(double));
    U = malloc((l_grid+2)*(l_grid+2)*3*sizeof(double));
    sendArray = malloc(l_grid*l_grid*sizeof(double));

    /* initiate the send/recv buffers */
    double sendVectorToRight[l_grid*3];   // @TODO: maybe we do not need the diagonal corners, so maybe we need to l_grid+2 for corners?, i.e. change the @TODO:corners
    double sendVectorToLeft[l_grid*3];
    double sendVectorToTop[l_grid*3];
    double sendVectorToBottom[l_grid*3];
    double recvVectorFromRight[l_grid*3];
    double recvVectorFromLeft[l_grid*3];
    double recvVectorFromTop[l_grid*3];
    double recvVectorFromBottom[l_grid*3];

    /* Initiate the shockwave (water height */    // seed a shockwave by raising the height of water at bottom right corner of the grid
    //printf("mex: %d, mey: %d rank: %d\n", mex, mey, rank);  // @andy:debug:2210
    //printf("x: %d, y: %d\n", physx(0, mex, l_grid), physy(0, mey, l_grid)); // @andy:debug:2210
    //printf("l_grid: %d\n", l_grid); // @andy:debug:2210

    for (int x = 1; x < l_grid+1; ++x) // x = 0 and x = l_grid+2 are unallocated "ghost layers", we skip these points by starting at x=1 until x=grid+1
    {
        for (int y = 1; y < l_grid+1; ++y)
        {
            // if (            // if (x<4 && y<4){
            //     printf("ANDY: y: %i, x: %i\n", y, x);    
            // }x<4 && y<4){
            //     printf("ANDY: y: %i, x: %i\n", y, x);    
            // }
            h[x*(l_grid+2)  + y] = 0.1;
            u[x*(l_grid+2)  + y] = 0.0;
            v[x*(l_grid+2)  + y] = 0.0;
            U[(x*(l_grid+2) + y)*3]     = h[x*(l_grid+2) + y];
            U[(x*(l_grid+2) + y)*3 + 1] = u[x*(l_grid+2) + y] * h[x*(l_grid+2) + y];
            U[(x*(l_grid+2) + y)*3 + 2] = v[x*(l_grid+2) + y] * h[x*(l_grid+2) + y];
            /* initialise the shock wave by lifting up parts of the water */
            if (   physx(x, mex, (l_grid)) < (xmax+1) 
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
        for (int y = 1; y < l_grid+2; ++y)
        {
            for (int i = 0; i < 3; ++i)
            {
                F[ (x*(l_grid + 2) + y)*3 + i ] = 0.0;
            }
        }
    }
    for (int x = 1; x < l_grid+2; ++x)
    {
        for (int y = 1; y < l_grid+2; ++y)
        {
            for (int i = 0; i < 3; ++i)
            {
                G[ (x*(l_grid + 2) + y)*3 + i ] = 0.0;
            }
        }
    }
    // calculate displacement before for loop dicrease computational effort, to properly organise the tiles into global grid
    int sendcounts[npx*npx];
    int displs[npx*npx];

    if (rank == 0) {
        for (int i=0; i<npx*npx; i++) sendcounts[i] = 1;
        int disp = 0;
        for (int i=0; i<npx; i++) {
            for (int j=0; j<npx; j++) {
                displs[i*npx+j] = disp;
                disp += 1;
            }
            disp += ((n_grid/npx-1))*npx;
        }
    }

    for (int i = 0; i < totalNumberofTimeStep; ++i)
    {
        // plot before sending
        if (rank == 1 && i == 10)
        {
            y = 1;
            for (int x = 1; x < l_grid+1; ++x)
            {
                printf("%lf ", U[ ((x)*(l_grid+2) + (y))*3 ]);
            }
            printf("\n");
        }

        // plot before sending
        if (rank == 0 && i == 10)
        {
            y = l_grid+1;
            for (int x = 1; x < l_grid+1; ++x)
            {
                printf("%lf ", U[ ((x)*(l_grid+2) + (y))*3 ]);
            }
            printf("\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
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

        /* Ghost layers "halo": transfer data into halo prior to computation */
        // @pseudo: copy stuff into sendbuffer
        // @pseudo:1 ^ data in the horizontal direction first
        // @pseudo:2 ^ verticle

        if(mex != (npx-1))    // mex=(npx-1) are the processors on the rightermost boundary (i.e. cannot send to the right, or recv from right)            
        {
            x = l_grid;
            for (int y = 0; y < l_grid; ++y)
            {
                debug = 1;
                if (debug==1){ printf("\trank: %i\t,time: %i\t,debug: %i\t, l_grid: %i\t, y: %i\n",rank,i,debug,l_grid,y); };  // @ANDY:debug:1 @debug1
                sendVectorToRight[y*3]     = U[ ((x)*(l_grid+2) + (y+1))*3 ];  // @TODO:corners
                sendVectorToRight[y*3 + 1] = U[ ((x)*(l_grid+2) + (y+1))*3 + 1 ]; 
                sendVectorToRight[y*3 + 2] = U[ ((x)*(l_grid+2) + (y+1))*3 + 2 ]; 
            }

        }

        if(mex != 0)    
        {     
            x = 1;
            for (int y = 0; y < l_grid; ++y)
            {
                sendVectorToLeft[y*3]     = U[ ((x)*(l_grid+2) + (y+1))*3 ];  
                sendVectorToLeft[y*3 + 1] = U[ ((x)*(l_grid+2) + (y+1))*3 + 1 ];   
                sendVectorToLeft[y*3 + 2] = U[ ((x)*(l_grid+2) + (y+1))*3 + 2];   
            }            
        }

        if(mey != (npx-1))    
        {     
            y = l_grid;
            for (int x = 0; x < l_grid; ++x)
            {
                sendVectorToTop[x*3]     = U[ ((x+1)*(l_grid+2) + (y))*3 ]; 
                sendVectorToTop[x*3 + 1] = U[ ((x+1)*(l_grid+2) + (y))*3 + 1 ]; 
                sendVectorToTop[x*3 + 2] = U[ ((x+1)*(l_grid+2) + (y))*3 + 2 ]; 
            }
        }

        if(mey != 0)    
        {
            y = 1;
            for (int x = 0; x < l_grid; ++x)
            {
                sendVectorToBottom[x*3]     = U[ ((x+1)*(l_grid+2) + (y))*3 ];    
                sendVectorToBottom[x*3 + 1] = U[ ((x+1)*(l_grid+2) + (y))*3 + 1];    
                sendVectorToBottom[x*3 + 2] = U[ ((x+1)*(l_grid+2) + (y))*3 + 2];    
            }            
        }

        // @pseudo: send sendbuffer into halo        
        // @pseudo:1 send data in the horizontal direction first
        if(mex != (npx-1))    // mex=(npx-1) are the processors on the rightermost boundary (i.e. cannot send to the right, or recv from right)            
        {     
            // MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
            MPI_Send( &(sendVectorToRight[0]), l_grid*3, MPI_DOUBLE, (mey+(mex+1)*npx), 1, MPI_COMM_WORLD);
            //                                                      ^ this gives you the target rank, given the mex and mey
            //                                                                        ^ tag: e.g. u/ to make a correspondence between a send and a recv, a send passes data only to a recv with the same tag 
            MPI_Recv( &(recvVectorFromRight[0]), l_grid*3, MPI_DOUBLE, (mey+(mex+1)*npx), 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //                                                        ^ this gives you the source rank, given the mex and mey, i.e. which rank to recv data from
        }

        if(mex != 0)    
        {     
            MPI_Send( &(sendVectorToLeft[0]), l_grid*3, MPI_DOUBLE, (mey+(mex-1)*npx), 2, MPI_COMM_WORLD);
            MPI_Recv( &(recvVectorFromLeft[0]), l_grid*3, MPI_DOUBLE, (mey+(mex-1)*npx), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // @pseudo:2 ^ verticle

        if(mey != (npx-1))    
        {
            MPI_Send( &(sendVectorToTop[0]), l_grid*3, MPI_DOUBLE, (mey+1+mex*npx), 3, MPI_COMM_WORLD);
            MPI_Recv( &(recvVectorFromTop[0]), l_grid*3, MPI_DOUBLE, (mey+1+mex*npx), 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        
        }

        if(mey != 0)
        {
            MPI_Send( &(sendVectorToBottom[0]), l_grid*3, MPI_DOUBLE, (mey-1+(mex)*npx), 4, MPI_COMM_WORLD);
            MPI_Recv( &(recvVectorFromBottom[0]), l_grid*3, MPI_DOUBLE, (mey-1+(mex)*npx), 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // @pseudo: copy stuff from recvbuffers into the U grid
        // @pseudo:1 ^ data in the horizontal direction first
        // @pseudo:2 ^ verticle
        if(mex != (npx-1))    // mex=(npx-1) are the processors on the rightermost boundary (i.e. cannot send to the right, or recv from right)            
        {     
            x = l_grid+1;
            for (int y = 0; y < l_grid; ++y)
            {
                U[ ((x)*(l_grid+2) + (y+1))*3 ]     = recvVectorFromRight[y*3];  // @TODO:corners
                U[ ((x)*(l_grid+2) + (y+1))*3 + 1 ] = recvVectorFromRight[y*3 + 1]; 
                U[ ((x)*(l_grid+2) + (y+1))*3 + 2 ] = recvVectorFromRight[y*3 + 2]; 
            }
        }

        if(mex != 0)    
        {     
            x = 0;
            for (int y = 0; y < l_grid; ++y)
            {
                U[ ((x)*(l_grid+2) + (y+1))*3 ]     = recvVectorFromLeft[y*3];  // @TODO:corners
                U[ ((x)*(l_grid+2) + (y+1))*3 + 1 ] = recvVectorFromLeft[y*3 + 1]; 
                U[ ((x)*(l_grid+2) + (y+1))*3 + 2 ] = recvVectorFromLeft[y*3 + 2]; 
            }            
        }

        if(mey != (npx-1))    
        {     
            y = l_grid+1;
            for (int x = 0; x < l_grid; ++x)
            {
                U[ ((x+1)*(l_grid+2) + (y))*3 ]     = recvVectorFromTop[x*3];  // @TODO:corners
                U[ ((x+1)*(l_grid+2) + (y))*3 + 1 ] = recvVectorFromTop[x*3 + 1]; 
                U[ ((x+1)*(l_grid+2) + (y))*3 + 2 ] = recvVectorFromTop[x*3 + 2]; 
            }
        }

        if(mey != 0)    
        {     
            y = 0;
            for (int x = 0; x < l_grid; ++x)
            {
                U[ ((x+1)*(l_grid+2) + (y))*3 ]     = recvVectorFromBottom[x*3];  // @TODO:corners
                U[ ((x+1)*(l_grid+2) + (y))*3 + 1 ] = recvVectorFromBottom[x*3 + 1]; 
                U[ ((x+1)*(l_grid+2) + (y))*3 + 2 ] = recvVectorFromBottom[x*3 + 2]; 
            }            
        }

        // plot after sending
        if (rank == 0 && i == 10)
        {
            y = l_grid+1;
            for (int x = 1; x < l_grid+1; ++x)
            {
                printf("%lf ", U[ ((x)*(l_grid+2) + (y))*3 ]);
            }
            printf("\n");
        }
        /* compute fluxes*/
        computeFlux(U, F, G, l_grid, &amax, mex, mey, npx);

        /* updating the fluxes*/
        updateFlux(U, F, G, l_grid, dt_dx);

        //printf("Time Step = %d, amax = %lf \n", i, amax);
        //printf("Time Step = %d, Courant Number = %lf \n", i, amax * dt_dx* 2 );
        /* gather all information from all processors to main */
        int counter = 0;
       
        if (i % plottingStep == 0)
        {
            /* initialise send buffer */
            for (int x = 0; x < l_grid; ++x) // so x = 0 and x = l_grid+2 are unallocated "ghost layers" 
            {
                for (int y = 0; y < l_grid; ++y)
                {
                    sendArray[ (x*l_grid + y) ] = U[ ((x+1)*(l_grid + 2) + (y+1))*3 ];
                }
            }
            // if (rank == 1)
            // {
            //     // write_vtkFile(szProblem, i, l_grid+2, l_grid+2, l_grid+2, cellsize, cellsize, h);
            //     write_vtkFile(szProblem, i, l_grid, l_grid, l_grid, cellsize, cellsize, sendArray);
            // }
            /* gather all tiles in root branch */
            MPI_Gatherv(&(sendArray[0]), l_grid*l_grid,  MPI_DOUBLE, U_global, sendcounts, displs, subarrayType, 0, MPI_COMM_WORLD);
            /* write vtk file*/
            if (rank == 0)
            {
                write_vtkFile(szProblem, i, length, n_grid, n_grid, cellsize, cellsize, U_global);
            }
        }
    }
    // free mpi data type
    MPI_Type_free(&subarrayType);

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
