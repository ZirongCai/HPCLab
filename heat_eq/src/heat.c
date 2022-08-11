#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "input.h"
#include "heat.h"
#include "timing.h"
//#include "misc.h"
#include "mmintrin.h"

#include "omp.h"
//#include <mpi.h>

//#include <papi.h>

// Comment out, if you don't want to compute the norm over all ranks
#define COMPUTE_NORM

#define NUM_EXPERIMENTS (2)

double* time;

void usage(char* s)
{
    fprintf(stderr, "Usage: %s <input file> [<N> <M>] [result file]\n\n", s);
}

double meand(double* a, int n)
{
    double sum = 0;
    for (int i = 0; i < n; ++i)
        sum += a[i];
    return sum / n;
}

int main(int argc, char* argv[])
{

    // algorithmic parameters
    algoparam_t param;

    // timing

    double residual;

    // set the visualization resolution
    param.visres = 100;

    // check arguments
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }
    // if number of arguments is not 2, 3, 4 or 5 error!
    // - Two if only inputfile is given
    // - 3 if input and outputfile is given
    // - 4 if input file and topology dims are given
    // - 5 if everything is given
    if (argc > 5) {
        usage(argv[0]);
        return 1;
    }

    // check input file
    FILE* infile = NULL;
    if (!(infile = fopen(argv[1], "r"))) {
        fprintf(stderr, "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);

        usage(argv[0]);
        return 1;
    }

    // check n and m only possible if there are 4 or 5 arguments
    int n = -1;
    int m = -1;
    if (argc == 4 || argc == 5) {
        n = atoi(argv[2]);
        m = atoi(argv[3]);
    }

    // check result file
    char* resfilename = (argc == 3) ? argv[2] : (argc == 5 ? argv[4] : "heat.ppm");

    FILE* resfile = NULL;
    if (!(resfile = fopen(resfilename, "w"))) {
        fprintf(stderr, "\nError: Cannot open \"%s\" for writing.\n\n", resfilename);

        usage(argv[0]);
        return 1;
    }

    // check input
    if (!read_input(infile, &param)) {
        fprintf(stderr, "\nError: Error parsing input file.\n\n");

        usage(argv[0]);
        return 1;
    }

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm;
    int dim[2], period[2], reorder;
    int coord[2], id;

    // Set rank
    param.rank = rank;

    // If n and m is set, override param
    if (n > 0 && m > 0) {
        param.dimx = n;
        param.dimy = m;

        if (rank == ROOT_RANK) {
            fprintf(stdout, "Using topology info from command line\n");
        }
    }

    if (rank == ROOT_RANK) {
        if (param.dimx < 1 || param.dimx < 1) {
            fprintf(stdout, "Can't work with topologies with sizes smaller than 1, given %d, %d\n", param.dimx,
                    param.dimy);
            return 1;
        }
        if (n == 0 || m == 0) {
            // Something went wrong parsing, or maybe wrong numbers
            fprintf(stderr, "Topology dims from command line can be considered, falling back to input file\n");
        }
    }

    // Print out input dimx dimy
    if (rank == ROOT_RANK) {
#ifndef BLOCKING_COMMUNICATION 
        fprintf(stdout, "Using non-blocking communication\n");
#else
        fprintf(stdout, "Using blocking communication\n");
#endif
    }

    // Check if we have enougth processes to split the wanted topo
    if (size != param.dimx * param.dimy) {
        fprintf(stdout, "Number of MPI processes (%d) not sufficient to created wanted topology (%d, %d)", size,
                param.dimx, param.dimy);
        MPI_Finalize();
        return 0;
    }

    // if we size > param.dimx * param.dimy, we can not rely on size of our topology calculation
    int act_size = param.dimx * param.dimy;

    split_topography(&param, &comm);

    if (rank == ROOT_RANK)
        time = (double*) calloc(sizeof(double),
                                (int) (param.max_res - param.initial_res + param.res_step_size) / param.res_step_size);

#if 0
    if (param.topo_x == 0) {
        fprintf(stdout, "[rank %d] I'm top row, so top should be negative -> %d\n", rank, param.topo_top);
    }

    if (param.topo_x == param.dimx - 1) {
        fprintf(stdout, "[rank %d] I'm bottom row, so bottom should be negative -> %d\n", rank, param.topo_down);
    }

    if (param.topo_y == 0) {
        fprintf(stdout, "[rank %d] I'm left column, so left should be negative -> %d\n", rank, param.topo_left);
    }

    if (param.topo_y == param.dimy - 1) {
        fprintf(stdout, "[rank %d] I'm right column, so right should be negative -> %d\n", rank, param.topo_right);
    }
#endif

    // Set rank, from here we need it in the struct
    param.rank = rank;

    if (rank == ROOT_RANK) {
        print_params(&param);
#if 0
        // Print some info to ensure we keep the results correct
        fprintf(stdout, "\nAll the values are taken from assignment 3, just as a reference if we mess up something\n");
        fprintf(stdout, "Res  500, residual 0.082628, norm 26.11\n");
        fprintf(stdout, "Res 1000, residual 0.165178, norm 36.91\n");
        fprintf(stdout, "Res 1500, residual 0.247708, norm 45.20\n");
        fprintf(stdout, "Res 2000, residual 0.330234, norm 52.19\n");
        fprintf(stdout, "Res 2500, residual 0.412757, norm 58.34\n");
        fprintf(stdout, "Res 3000, residual 0.495280, norm 63.91\n");
        fprintf(stdout, "Res 3500, residual 0.577802, norm 69.03\n");
        fprintf(stdout, "Res 4000, residual 0.660324, norm 73.80\n");
        fprintf(stdout, "Res 4500, residual 0.742846, norm 78.27\n\n\n");
#endif

        fprintf(stdout, "Domain distribution %d %d \n", param.dimx, param.dimy);
#ifdef COMPUTE_NORM
        fprintf(stdout, "| %8s | %8s | %8s | %10s | %6s |\n", "Res", "Time", "Residual", "GFlops", "Norm");
        fprintf(stdout, "|----------|----------|----------|------------|--------|\n");
#else
        fprintf(stdout, "| %8s | %8s | %8s | %10s |\n", "Res", "Time", "Residual", "GFlops");
        fprintf(stdout, "|----------|----------|----------|------------|\n");
#endif
    }

    int exp_number = 0;

    double* ufinal = NULL;

    for (param.act_res = param.initial_res; param.act_res <= param.max_res;
         param.act_res = param.act_res + param.res_step_size) {

        if (!init_topology(&param)) {
            fprintf(stderr, "Couldn't split into desired topology\n");

            usage(argv[0]);
            break;
        }

        const int np = param.act_res + 2;

        // nx and ny plus ghost cells
        const int nx = param.topo_nx + 2;
        const int ny = param.topo_ny + 2;

        // Some variables to store and print everything in
        // We only store the last residual and norm, hopefully we have deterministic behaviour
        double last_residual = 0.0;
        double last_norm = 0.0;

        // Arrays to store time and bandwidth for averaging it
        double atime[NUM_EXPERIMENTS];
        double abandwidth[NUM_EXPERIMENTS];

        for (int exp = 0; exp < NUM_EXPERIMENTS; ++exp) {
            if (!initialize(&param)) {
                fprintf(stderr, "Error in Jacobi initialization.\n\n");

                usage(argv[0]);
            }

            // Copy over NUMA optimized
#pragma omp parallel for schedule(guided)
            for (int i = 0; i < ny; i++) {
                for (int j = 0; j < nx; j++) {
                    param.uhelp[i * nx + j] = param.u[i * nx + j];
                }
            }

            // david.frank TODO: do we need that?
            MPI_Barrier(MPI_COMM_WORLD);

            // starting time
            if (rank == ROOT_RANK)
                time[exp_number] = wtime();

            residual = 999999999;

            // All the iteration loops are no in their seperate function
            double tmp_residual = jacobi_iters(&param);

            // We do the residual reduce outside of the loop and only at the end.
            // We would need to do this earler inside of the loop, if we had any stopping criteria, but
            // as we don't have it, we should are fine only reducing once we stop looping
            // Note that the results are the same, but we don't have to sync on it every iteration, if we don't need it
            MPI_Reduce(&tmp_residual, &residual, 1, MPI_DOUBLE, MPI_SUM, ROOT_RANK, MPI_COMM_WORLD);

            double t = 0.0;
            if (rank == ROOT_RANK) {
                time[exp_number] = wtime() - time[exp_number];
                t = time[exp_number];
            }

#ifdef COMPUTE_NORM
            double tmp_norm = 0.0;
            for (int i = 1; i < ny - 1; ++i) {
                for (int j = 1; j < nx - 1; j++) {
                    int index = i * nx + j;
                    tmp_norm += param.u[index] * param.u[index];
                }
            }

            double norm = 0;
            MPI_Reduce(&tmp_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            last_norm = sqrt(norm);
#endif
            // In all iterations but the last get ride of u, uhelp and uvis
            if (exp < NUM_EXPERIMENTS - 1) {
                finalize(&param);
            } else {
                if (ufinal) {
                    free(ufinal);
                    ufinal = NULL;
                }
                ufinal = param.u;
                param.u = NULL;

                if (param.uhelp) {
                    free(param.uhelp);
                    param.uhelp = NULL;
                }
            }

            last_residual = residual;
            /** atime[exp] = time[exp_number]; */
            /** abandwidth[exp] = (double) param.maxiter * (np - 2) * (np - 2) * 7 / time[exp_number] / 1000000000; */
            double gflops = (double) param.maxiter * (np - 2) * (np - 2) * 7 / t / 1000000000;
            abandwidth[exp] = gflops;
            atime[exp] = t;
        }

        // Calculate mean
        const double avg_time = meand(atime, NUM_EXPERIMENTS);
        const double avg_bandwidth = meand(abandwidth, NUM_EXPERIMENTS);

        if (rank == ROOT_RANK) {

#ifdef COMPUTE_NORM
            fprintf(stdout, "| %8d | %8f | %8f | %10.3f | %6.2f |\n", param.act_res, avg_time, last_residual,
                    avg_bandwidth, last_norm);
#else
            fprintf(stdout, "| %8d | %8f | %8f | %10.3f |\n", param.act_res, atime[exp], last_residual,
                    abandwidth[exp]);
#endif

            /** fprintf(stdout, "| %4d | %10d | %2d x %2d | (%2d, %2d) | %4d x %4d | (%4d,%4d) | (%4d,%4d) |\n",
             * rank,
             */
            /**         param.act_res, param.dimx, param.dimy, param.topo_x, param.topo_y, param.topo_nx,
             * param.topo_ny,
             */
            /**         param.glob_startx, param.glob_starty, param.glob_endx, param.glob_endy); */

            // Debug print the thing once
            if (act_size == 1 && param.act_res <= 20) {
                // Printf it once
                // fprintf(stdout, "");
                for (int i = 1; i < ny - 1; ++i) {
                    for (int j = 1; j < nx - 2; ++j) {
                        int index = i * nx + j;
                        fprintf(stdout, "%8.5f, ", param.u[index]);
                    }
                    fprintf(stdout, "%8.5f\n", param.u[i * nx + nx - 2]);
                }
            }
        }

        exp_number++;
    }

    // Maybe get this to work at some point again, but reeeeealy don't care to much about it now :D
    /** param.act_res = param.act_res - param.res_step_size; */

    //write_image(resfile, ufinal, param.topo_nx, param.topo_ny);
    // ufinal has the content of the last iteration of the loop 
//    coarsen(ufinal, param.topo_nx, param.topo_ny, param.uvis, param.visres + 2, param.visres + 2);
    //write_image(resfile, param.uvis, param.visres+2, param.visres+2);
    coarsen(ufinal, &param);
//    add_dummy_data(&param);

    unsigned int dummy_resx = param.visres * param.dimx;
    unsigned int dummy_resy = param.visres * param.dimy;
    double* global_data;
    if (rank == 0)
        global_data = (double *) malloc(dummy_resx * dummy_resy * sizeof(double));
    
    all_gatherv(&param, &comm, size, rank, global_data);
   
     
    // Don't forget to free ufinal as well
    free(ufinal);
    free(param.uvis);
   
//    unsigned int global_res[2];
//    MPI_Reduce(&param.lresx, &global_res[0], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&param.lresy, &global_res[1], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
            //write_image(resfile, global_data, global_res[0], global_res[1], dummy_res);
            //write_image(resfile, global_data, dummy_res, dummy_res, dummy_res);
            write_image(resfile, global_data, dummy_resx, dummy_resy);
            free(global_data);
    }
    if (rank==0)
        free(time);
    free(param.heatsrcs);
    MPI_Finalize();
    return 0;
}

