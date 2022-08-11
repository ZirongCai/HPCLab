#pragma once

#include "heat.h"

/// Blocking communication between MPI processes
int blocking_communication(algoparam_t* param, MPI_Datatype column_type, MPI_Request* reqs);

/// Non blocking communication between MPI processes
int nonblocking_communication(algoparam_t* param, MPI_Datatype column_type, MPI_Request* reqs);
