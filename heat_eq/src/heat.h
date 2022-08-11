/*
 * heat.h
 *
 * Global definitions for the iterative solver
 */

#ifndef JACOBI_H_INCLUDED
#define JACOBI_H_INCLUDED

#include <stdio.h>
#include <mpi.h>

#define MY_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MY_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define ROOT_RANK (0)

// #define BLOCKING_COMMUNICATION 1
//#define RUN_HYBRID
#define COMPUTE_WHILE_WAITING

// configuration

typedef struct {
    float posx;
    float posy;
    float range;
    float temp;
} heatsrc_t;

typedef struct {
    unsigned maxiter; // maximum number of iterations
    unsigned act_res;
    unsigned max_res; // spatial resolution
    unsigned initial_res;
    unsigned res_step_size;
    unsigned visres; // visualization resolution

    double *u, *uhelp;
    double* uvis;
    int dimx, dimy;

    /// Size which the current block works on (no boundary cells considered)
    int topo_nx;
    int topo_ny;

    /// Store ranks x and y coord in topology
    int topo_x, topo_y;

    /// Store ranks top and bottom neighbour in topology
    int topo_top, topo_down;

    /// Store ranks left and right neighbour in topology
    int topo_left, topo_right;

    /// Coords for global domain [start, end[ (end is exclusiv)
    int glob_startx, glob_starty;
    int glob_endx, glob_endy;

    /// Rank of current process
    int rank;

    /// Local resolution
    unsigned int lresx, lresy;

    /// average local nx, ny (for processes that are not located at last the row or column in vitual topo)
    int avrg_nx, avrg_ny;

    unsigned numsrcs; // number of heat sources
    heatsrc_t* heatsrcs;
} algoparam_t;

// function declarations

// misc.c
int init_topology(algoparam_t* param);
int initialize(algoparam_t* param);
int finalize(algoparam_t* param);
//void write_image(FILE* f, double* u, unsigned sizex, unsigned sizey, unsigned int dummy_size);
void write_image(FILE* f, double* u, unsigned sizex, unsigned sizey);
//int coarsen(double* uold, unsigned oldx, unsigned oldy, double* unew, unsigned newx, unsigned newy);
int coarsen(double* uold,algoparam_t* param);
// Gauss-Seidel: relax_gauss.c
double residual_gauss(double* u, double* utmp, unsigned sizex, unsigned sizey);
void relax_gauss(double* u, unsigned sizex, unsigned sizey);

// Jacobi: relax_jacobi.c
double residual_jacobi(double* u, unsigned sizex, unsigned sizey);
double relax_jacobi(double** u, double** utmp, unsigned sizex, unsigned sizey);
double relax_jacobi_hybrid(double** u, double** utmp, unsigned sizex, unsigned sizey);
void split_topography(algoparam_t* param, MPI_Comm* comm);
double jacobi_iters(algoparam_t* param);
double outer_jacobi(double* u, double* utmp1, unsigned sizex, unsigned sizey, unsigned borderx, unsigned bordery);
double inner_jacobi(double* u, double* utmp, unsigned sizex, unsigned sizey, unsigned offsetx, unsigned offsety);
void add_dummy_data(algoparam_t* param);
void all_gatherv(algoparam_t* param, MPI_Comm* comm, int size, int rank, double* gblock_ptr);
#endif // JACOBI_H_INCLUDED
