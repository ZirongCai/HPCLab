/*
 * misc.c
 *
 * Helper functions for
 * - initialization
 * - finalization,
 * - writing out a picture
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "heat.h"

int init_topology(algoparam_t* param)
{
    // Figure out the size for each block
    const int dimx = param->dimx;
    const int dimy = param->dimy;

    // total number of points (including border)
    const int tmp_np = param->act_res;

    // Split each domain into dimx and dimy equal parts, don't consider boundary cells yet
    const int tmp_nx = ((tmp_np + dimx - 1) / dimx);
    const int tmp_ny = ((tmp_np + dimy - 1) / dimy);

    // find start and end coord in domain for each block in topology
    // start at 1, add nx * topo_x then subtract 1 for ghost cell at the left/top
    const int startx = 1 + (param->topo_x * tmp_nx) - 1;
    const int starty = 1 + (param->topo_y * tmp_ny) - 1;

    // Add two at the end, as wee start at 1, plus 1 ghost cell
    const int endx = MY_MIN((param->topo_x + 1) * tmp_nx + 2, tmp_np + 2);
    const int endy = MY_MIN((param->topo_y + 1) * tmp_ny + 2, tmp_np + 2);

    // Store average nx ny for later uses - Phuong
    param->avrg_nx = ceil(tmp_nx);
    param->avrg_ny = ceil(tmp_ny);

    param->glob_startx = startx;
    param->glob_starty = starty;

    param->glob_endx = endx;
    param->glob_endy = endy;

    param->topo_nx = endx - startx - 2;
    param->topo_ny = endy - starty - 2;

    return 1;
}

/*
 * Initialize the iterative solver
 * - allocate memory for matrices
 * - set boundary conditions according to configuration
 */
int initialize(algoparam_t* param)
{
    double dist;

    init_topology(param);

    const int nx = param->topo_nx + 2;
    const int ny = param->topo_ny + 2;

    const int startx = param->glob_startx;
    const int starty = param->glob_starty;

    const int endx = param->glob_endx;
    const int endy = param->glob_endy;

    /** if (rank == 0) */
    /**     fprintf(stdout, "[rank %d] dimx = %d, dimy = %d\n", param->rank, dimx, dimy); */
    /** if (rank == 0) */
    /** fprintf(stdout, "[rank %d] nx = %d, ny = %d, tmp_nx %d, tmp_ny %d --- starty %d, endy %d\n\n\n", param->rank,
     * nx, ny, */
    /**         tmp_nx, tmp_ny, starty, endy); */

    /** fprintf(stdout, "[rank %d] startx = %d, starty = %d\n", param->rank, startx, starty); */
    /** fprintf(stdout, "[rank %d] endx = %d, endy = %d\n", param->rank, endx, endy); */

    //
    // allocate memory
    //
    // Assume that we store row major, so access row i and col j => i * nx + j
    (param->u) = (double*) malloc(sizeof(double) * nx * ny);
    (param->uhelp) = (double*) malloc(sizeof(double) * nx * ny);
    (param->uvis) = (double*) calloc(sizeof(double), (param->visres + 2) * (param->visres + 2));

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            param->u[i * nx + j] = 0;
            param->uhelp[i * nx + j] = 0;
        }
    }

    if (!(param->u) || !(param->uhelp) || !(param->uvis)) {
        fprintf(stderr, "Error: Cannot allocate memory\n");
        return 0;
    }

    const int np = param->act_res + 2;

    // For each heat source
    for (int i = 0; i < param->numsrcs; i++) {
        /* top row */
        if (param->topo_y == 0) {
            // Only work on the part of x of your block, keep a local index running to access array (lj)
            for (int j = startx, lj = 0; j < endx; j++, lj++) {
                // Calculate distance based on global coordinates
                dist = sqrt(pow((double) j / (double) (np - 1) - param->heatsrcs[i].posx, 2)
                        + pow(param->heatsrcs[i].posy, 2));

                /** dist = sqrt(pow(param->heatsrcs[i].posx, 2) */
                /**             + pow((double) j / (double) (np - 1) - param->heatsrcs[i].posy, 2)); */

                // Set value if we're in range of heat source
                if (dist <= param->heatsrcs[i].range) {
                    (param->u)[lj] +=
                        (param->heatsrcs[i].range - dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
                }
            }
        }

        /* bottom row, indicated by topo_y == n-1 */
        if (param->topo_y == param->dimy - 1) {
            // Iterate over your part of the row
            for (int j = startx, lj = 0; j < endx; j++, lj++) {
                // Calculate distance based on global coordinates
                dist = sqrt(pow((double) j / (double) (np - 1) - param->heatsrcs[i].posx, 2)
                        + pow(1 - param->heatsrcs[i].posy, 2));

                if (dist <= param->heatsrcs[i].range) {
                    // Access has to bo local again
                    // (ny - 1) * nx is the index for the last row
                    (param->u)[(ny - 1) * nx + lj] +=
                        (param->heatsrcs[i].range - dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
                }
            }
        }

        // we don't handle the first and last row, they are handled above
        const int start = starty;
        const int end = endy;

        /* leftmost column */
        if (param->topo_x == 0) {
            for (int j = start, lj = 0; j < end; j++, lj++) {
                dist = sqrt(pow(param->heatsrcs[i].posx, 2)
                        + pow((double) j / (double) (np - 1) - param->heatsrcs[i].posy, 2));

                if (dist <= param->heatsrcs[i].range) {
                    (param->u)[lj * nx] +=
                        (param->heatsrcs[i].range - dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
                }
            }
        }

        /* rightmost column */
        if (param->topo_x == param->dimx - 1) {
            for (int j = start, lj = 0; j < end; j++, lj++) {
                dist = sqrt(pow(1 - param->heatsrcs[i].posx, 2)
                        + pow((double) j / (double) (np - 1) - param->heatsrcs[i].posy, 2));

                if (dist <= param->heatsrcs[i].range) {
                    (param->u)[lj * nx + (nx - 1)] +=
                        (param->heatsrcs[i].range - dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
                }
            }
        }
    }

    return 1;
}

/*
 * free used memory
 */
int finalize(algoparam_t* param)
{
    if (param->u) {
        free(param->u);
        param->u = 0;
    }

    if (param->uhelp) {
        free(param->uhelp);
        param->uhelp = 0;
    }

    if (param->uvis) {
        free(param->uvis);
        param->uvis = 0;
    }
    return 1;
}

void write_image(FILE* f, double* u, unsigned sizex, unsigned sizey)
{
    // RGB table
    unsigned char r[1024], g[1024], b[1024];
    int i, j, k;

    double min, max;

    j = 1023;

    // prepare RGB table
    for (i = 0; i < 256; i++) {
        r[j] = 255;
        g[j] = i;
        b[j] = 0;
        j--;
    }
    for (i = 0; i < 256; i++) {
        r[j] = 255 - i;
        g[j] = 255;
        b[j] = 0;
        j--;
    }
    for (i = 0; i < 256; i++) {
        r[j] = 0;
        g[j] = 255;
        b[j] = i;
        j--;
    }
    for (i = 0; i < 256; i++) {
        r[j] = 0;
        g[j] = 255 - i;
        b[j] = 255;
        j--;
    }

    min = DBL_MAX;
    max = -DBL_MAX;

    // find minimum and maximum
    for (i = 0; i < sizey; i++) {
        for (j = 0; j < sizex; j++) {
            if (u[i * sizex + j] > max)
                max = u[i * sizex + j];
            if (u[i * sizex + j] < min)
                min = u[i * sizex + j];
        }
    }

    fprintf(f, "P3\n");
    fprintf(f, "%u %u\n", sizex, sizey);
    fprintf(f, "%u\n", 255);

    for (i = 0; i < sizey; i++) {
        for (j = 0; j < sizex; j++) {
            k = (int) (1024.0 * (u[i * sizex + j] - min) / (max - min));
            fprintf(f, "%d %d %d  ", r[k], g[k], b[k]);
        }
        fprintf(f, "\n");
    }
}

/*
//int coarsen(double* uold, unsigned oldx, unsigned oldy, double* unew, unsigned newx, unsigned newy)
int coarsen(double *uold, algoparam_t* param)
{
    int i, j, k, l, ii, jj;
    unsigned int newx = param->visres+2, newy=newx;
	// minus 2 = not including boundaries
    unsigned oldx = param->topo_nx-2, oldy=param->topo_ny-2;
    double* unew = param->uvis;

	double** sub_u = (double**) malloc(oldx* sizeof(double));
	    for (int i = 1; i < oldx; i++)
	        sub_u[i-1] = &(uold[i*(param->topo_ny)+1]);

	printf("Here 1 \n");
    // Recalculate resolution scale for processes that have less elements 
    if (oldx <  param->avrg_nx )
        newx = newx * oldx / param->avrg_nx;
    if (oldy <  param->avrg_ny )
        newy = newy * oldy / param->avrg_ny;

    int stopx = newx;
    int stopy = newy;
    // Store all local res for later uses
    param->lresx = ceil(stopx);
    param->lresy = ceil(stopy);




    float temp;
    float stepx = (float) oldx / (float) newx;
    float stepy = (float) oldy / (float) newy;

    if (oldx < newx) {
        stopx = oldx* oldx/param->avrg_nx;
        stepx = 1.0;
    }
    if (oldy < newy) {
        stopy = oldy* oldy/param->avrg_ny;
        stepy = 1.0;
    }

    // printf("oldx=%d, newx=%d\n",oldx,newx);
    // printf("oldy=%d, newy=%d\n",oldy,newy);
    // printf("rx=%f, ry=%f\n",stepx,stepy);
    // NOTE: this only takes the top-left corner,
    // and doesnt' do any real coarsening

    for (i = 0; i < stopx; i++) {
        ii = stepx * i;
        for (j = 0; j < stopy; j++) {
            jj = stepy * j;
            temp = 0;
            for (k = 0; k < stepx; k++) {
                for (l = 0; l < stepy; l++) {
                    if (ii + k < oldx && jj + l < oldy)
                        //temp += uold[(ii + k) * oldx + (jj + l)];
                        temp += sub_u[ii + k][jj + l];
                }
            }
            //unew[i * newx + j] = temp / (stepy * stepx);
            unew[j * (param->visres+2) + i] = temp / (stepy *stepx);
        }
    }
	free(sub_u);
    return 1;
}
*/
//int coarsen(double* uold, unsigned oldx, unsigned oldy, double* unew, unsigned newx, unsigned newy)
int coarsen(double *uold, algoparam_t* param)
{
    int i, j, k, l, ii, jj;
    
    unsigned int newx = param->visres+2, newy=newx;
	// minus 2 = not including boundaries
    unsigned oldx = param->topo_nx-2, oldy=param->topo_ny-2;
    double* unew = param->uvis;


    int stopx = newx;
    int stopy = newy;
    float temp;
    float stepx = (float) oldx / (float) newx;
    float stepy = (float) oldy / (float) newy;

    if (oldx < newx) {
        stopx = oldx;
        stepx = 1.0;
    }
    if (oldy < newy) {
        stopy = oldy;
        stepy = 1.0;
    }

    // printf("oldx=%d, newx=%d\n",oldx,newx);
    // printf("oldy=%d, newy=%d\n",oldy,newy);
    // printf("rx=%f, ry=%f\n",stepx,stepy);
    // NOTE: this only takes the top-left corner,
    // and doesnt' do any real coarsening

    for (i = 0; i < stopy; i++) {
        ii = stepy * i;
        for (j = 0; j < stopx; j++) {
            jj = stepx * j;
            temp = 0;
            for (k = 0; k < stepy; k++) {
                for (l = 0; l < stepx; l++) {
                    if (ii + k < oldx && jj + l < oldy)
                        temp += uold[(ii + k) * oldx + (jj + l)];
                }
            }
            unew[i * newx + j] = temp / (stepy * stepx);
        }
    }

    return 1;
}
