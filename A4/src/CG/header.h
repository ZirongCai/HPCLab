#include "mpi.h"

typedef struct {
//    unsigned maxiter; // maximum number of iterations
//    unsigned act_res;
//    unsigned max_res; // spatial resolution
//    unsigned initial_res;
//    unsigned res_step_size;
//    unsigned visres; // visualization resolution

//    double *u, *uhelp;
//    double* uvis;

    double * grid;
    double * b;
    std::size_t grid_points_1d;
    std::size_t actual_gridsize;
    double mesh_width;

    int dimx, dimy;

    /// Size which the current block works on (no boundary cells considered)
    std::size_t topo_nx, topo_ny; //with 2 ghosh layers
    std::size_t avrg_nx, avrg_ny; //without ghost layers

    /// Store ranks x and y coord in topology
    int topo_x, topo_y;

    /// Store ranks top and bottom neighbour in topology
    int topo_top, topo_bottom;

    /// Store ranks left and right neighbour in topology
    int topo_left, topo_right;

    /// Coords for global domain [start, end[ (end is exclusiv)
    std::size_t glob_startx, glob_starty;
    std::size_t glob_endx, glob_endy;

    /// Rank of current process
    int rank;

    /// Local resolution
    unsigned int lresx, lresy;

    int mpi_size;
    MPI_Comm comm_2d;

    /// average local nx, ny (for processes that are not located at last the row or column in vitual topo)    int avrg_nx, avrg_ny;

} struct_param;

/// store number of grid points in one dimension
//std::size_t grid_points_1d;

/// store begin timestep
//struct timeval begin;
/// store end timestep
//struct timeval end;

/**
 * initialize and start timer
 */
void timer_start(struct timeval &begin);

double timer_stop(struct timeval &begin, struct timeval &end);

void store_grid(struct_param *param, double *grid, std::string filename);

 double eval_init_func(double x, double y);

void init_grid(struct_param *param);

void init_b(struct_param* param);

//void g_copy(double* dest, double* src, std::size_t grid_points_1d);

double g_dot_product(struct_param* param, double* grid1, double* grid2);

void g_scale(struct_param *param, double* grid, double scalar);

void g_scale_add(struct_param *param, double* dest, double* src, double scalar);

void g_product_operator(struct_param *param, double* grid, double* result);
void solve(struct_param *param, std::size_t cg_max_iterations, double cg_eps);

