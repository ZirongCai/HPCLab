/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "blocks/SWE_Block.hh"
#include "SWE_BlockCudaAwareMPI.hh"
#include "SWE_BlockCUDA_kernels.hh"

#include "SWE_WavePropagationBlockCudaAwareMPI.hh"

#include "tools/help.hh"
#include "tools/Logger.hh"

#include <cassert>
#include <cstdlib>
#include <cmath>

//uncomment this to initialize in gpu
//#define INIT

using namespace std;

/*
 * helper function to read CUDA error codes
 * (implementation in swe.cu */
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "\nCuda error (%s): %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }
}

/*
 * helper function to read CUDA error codes
 * (implementation in swe.cu */
void tryCUDA(cudaError_t err, const char *msg)
{
    if( cudaSuccess != err)
    {
        fprintf(stderr, "\nCuda error (%s): %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }
}

SWE_Block* getCudaBlockInstance(float nx, float ny, float dx, float dy) {
  SWE_Block *block = new SWE_WavePropagationBlockCudaAwareMPI(nx, ny, dx, dy);
  return block;
}

/**
 * Constructor: allocate variables for simulation
 *
 * unknowns h,hu,hv,b are defined on grid indices [0,..,nx+1]*[0,..,ny+1]
 * -> computational domain is [1,..,nx]*[1,..,ny]
 * -> plus ghost cell layer
 *
 * flux terms are defined for edges with indices [0,..,nx]*[1,..,ny]
 * or [1,..,nx]*[0,..,ny] (for horizontal/vertical edges)
 * Flux term with index (i,j) is located on the edge between
 * cells with index (i,j) and (i+1,j) or (i,j+1)
 *
 * bathymetry source terms are defined for cells with indices [1,..,nx]*[1,..,ny]
 *
 *
 * @param i_cudaDevice ID of the CUDA-device, which should be used.
 */
SWE_BlockCudaAwareMPI::SWE_BlockCudaAwareMPI(
	int l_nx, int l_ny,
	float l_dx, float l_dy)
  : SWE_Block(l_nx, l_ny, l_dx, l_dy)
{
  if (nx % TILE_SIZE != 0) {
	cout << "WARNING: nx not a multiple of TILE_SIZE  -> will lead to crashes!"
	  << endl << flush;
  };
  if (ny % TILE_SIZE != 0) {
	cout << "WARNING: ny not a multiple of TILE_SIZE  -> will lead to crashes!"
	  << endl << flush;
  };
#ifdef PACKING_DATA
  // allocate consecutive memory for 2 columns with three unknowns each
  // (h, hu, hv, excluding b) for copy/ghost layer at bottom/top boundary
  int size = nx + 2;
  //bottomLayer = new float[6 * size];
  cudaMallocManaged((void**) &bottomLayer, 6 * size * sizeof(float));
  checkCUDAError("allocate Unified memory for bottom copy/ghost layer");
  bottomGhostLayer = new SWE_Block1D(bottomLayer, bottomLayer + size, bottomLayer + (2 * size), size);
  bottomCopyLayer = new SWE_Block1D(bottomLayer + (3 * size), bottomLayer + (4 * size), bottomLayer + (5 * size), size);

  // same for top boundary:
  //topLayer = new float[6 * size];
  cudaMallocManaged((void**) &topLayer, 6 * size * sizeof(float));
  checkCUDAError("allocate Unified memory for top copy/ghost layer");
  topGhostLayer = new SWE_Block1D(topLayer, topLayer + size, topLayer + (2 * size), size);
  topCopyLayer = new SWE_Block1D(topLayer + (3 * size), topLayer + (4 * size), topLayer + (5 * size), size);
#endif

}

void SWE_BlockCudaAwareMPI::initScenario(float _offsetX, float _offsetY,
                                 SWE_Scenario &i_scenario,
                                 const bool i_multipleBlocks)
{
  offsetX = _offsetX;
  offsetY = _offsetY;
#if defined(INIT)
    {
      float *hd;
      float *hud;
      float *hvd;
      float *bd;
      cudaMallocManaged(&hd, (nx + 2) * (ny + 2) * sizeof(float));
      cudaMallocManaged(&hud, (nx + 2) * (ny + 2) * sizeof(float));
      cudaMallocManaged(&hvd, (nx + 2) * (ny + 2) * sizeof(float));
      cudaMallocManaged(&hvd, (nx + 2) * (ny + 2) * sizeof(float));
      cudaMallocManaged(&bd, (nx + 2) * (ny + 2) * sizeof(float));

      for (int i = 1; i <= nx; i++)
      {
        for (int j = 1; j <= ny; j++)
        {
          float x = offsetX + (i - 0.5f) * dx;
          float y = offsetY + (j - 0.5f) * dy;
          hd[(i) * (ny + 2) + (j)] = i_scenario.getWaterHeight(x, y);
          hud[(i) * (ny + 2) + (j)] = i_scenario.getVeloc_u(x, y) * hd[(i) * (ny + 2) + (j)];
          hvd[(i) * (ny + 2) + (j)] = i_scenario.getVeloc_v(x, y) * hd[(i) * (ny + 2) + (j)];
        }
      }

      for (int i = 0; i <= nx + 1; i++)
      {
        for (int j = 0; j <= ny + 1; j++)
        {
          bd[i * (ny + 2) + j] = i_scenario.getBathymetry(offsetX + (i - 0.5f) * dx,
                                                          offsetY + (j - 0.5f) * dy);
        }
      }

      //initialize water height and discharge
      dim3 dimBlock(TILE_SIZE, TILE_SIZE);
      dim3 dimGrid(nx / TILE_SIZE, ny / TILE_SIZE);


      initScenario_cuda<<<dimGrid, dimBlock>>>(h.elemVector(), hu.elemVector(), hv.elemVector(), b.elemVector(), hd, hud, hvd, bd, nx, ny);

    }
#else
    {
      for (int i = 1; i <= nx; i++)
      {
        for (int j = 1; j <= ny; j++)
        {
          float x = offsetX + (i - 0.5f) * dx;
          float y = offsetY + (j - 0.5f) * dy;
          h[i][j] = i_scenario.getWaterHeight(x, y);
          hu[i][j] = i_scenario.getVeloc_u(x, y) * h[i][j];
          hv[i][j] = i_scenario.getVeloc_v(x, y) * h[i][j];
        }
      }

      // initialize bathymetry
      for (int i = 0; i <= nx + 1; i++)
      {
        for (int j = 0; j <= ny + 1; j++)
        {
          b[i][j] = i_scenario.getBathymetry(offsetX + (i - 0.5f) * dx,
                                             offsetY + (j - 0.5f) * dy);
        }
      }
    }
#endif


  // in the case of multiple blocks the calling routine takes care about proper boundary conditions.
  if (i_multipleBlocks == false)
  {
    // obtain boundary conditions for all four edges from scenario
    setBoundaryType(BND_LEFT, i_scenario.getBoundaryType(BND_LEFT));
    setBoundaryType(BND_RIGHT, i_scenario.getBoundaryType(BND_RIGHT));
    setBoundaryType(BND_BOTTOM, i_scenario.getBoundaryType(BND_BOTTOM));
    setBoundaryType(BND_TOP, i_scenario.getBoundaryType(BND_TOP));
  }

  // perform update after external write to variables
  synchAfterWrite();
}

/**
 * Destructor: de-allocate all variables
 */
SWE_BlockCudaAwareMPI::~SWE_BlockCudaAwareMPI() {
#ifdef PACKING_DATA
	cudaFree(topLayer);
	cudaFree(bottomLayer);
#endif
}

//==================================================================
// methods for simulation
//==================================================================

/**
 * set the values of all ghost cells depending on the specifed
 * boundary conditions
 */
void SWE_BlockCudaAwareMPI::setBoundaryConditions() {
#ifdef DBG
 cout << "Call kernel to compute h in ghost layer corner (for visualisation only) "
      << flush << endl;
#endif

#ifdef DBG
 cout << "Call kernel to compute left/right boundaries " << flush << endl;
#endif
//   synchWaterHeightAfterWrite();
//   synchDischargeAfterWrite();

  if (boundary[BND_LEFT] == PASSIVE || boundary[BND_LEFT] == CONNECT) {}
  else {
     dim3 dimBlock(1,TILE_SIZE);
     dim3 dimGrid(1,ny/TILE_SIZE);
     kernelLeftBoundary<<<dimGrid,dimBlock>>>(
        h.elemVector(),hu.elemVector(),hv.elemVector(),nx,ny,boundary[BND_LEFT]);
  };

  if (boundary[BND_RIGHT] == PASSIVE || boundary[BND_RIGHT] == CONNECT) {}
  else {
     dim3 dimBlock(1,TILE_SIZE);
     dim3 dimGrid(1,ny/TILE_SIZE);
     kernelRightBoundary<<<dimGrid,dimBlock>>>(
        h.elemVector(),hu.elemVector(),hv.elemVector(),nx,ny,boundary[BND_RIGHT]);
  };

#ifdef DBG
  cout << "Call kernel to compute bottom/top boundaries " << flush << endl;
#endif
  if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT) {}
  else
  {
      dim3 dimBlock(TILE_SIZE,1);
      dim3 dimGrid(nx/TILE_SIZE,1);
      kernelBottomBoundary<<<dimGrid,dimBlock>>>(
         h.elemVector(),hu.elemVector(),hv.elemVector(),nx,ny,boundary[BND_BOTTOM]);
  }
  if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT) {}
  else
  {
      dim3 dimBlock(TILE_SIZE,1);
      dim3 dimGrid(nx/TILE_SIZE,1);
      kernelTopBoundary<<<dimGrid,dimBlock>>>(
         h.elemVector(),hu.elemVector(),hv.elemVector(),nx,ny,boundary[BND_TOP]);
  }

  kernelHdBufferEdges<<<1,1>>>(h.elemVector(), nx, ny);
  kernelHdBufferEdges<<<1,1>>>(hu.elemVector(), nx, ny);
  kernelHdBufferEdges<<<1,1>>>(hv.elemVector(), nx, ny);

}

#ifdef PACKING_DATA
void SWE_BlockCudaAwareMPI::synchCopyLayerBeforeRead()
{
#ifdef DBG
  cout << "Packing data ..." << flush << endl;
#endif

  int size = 3 * (nx + 2);
  // bottom copy layer
  if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT)
  {
    dim3 dimBlock(TILE_SIZE, 1);
    dim3 dimGrid(nx / TILE_SIZE, 1);
    kernelBottomCopyLayer<<<dimGrid, dimBlock>>>(
        h.elemVector(), hu.elemVector(), hv.elemVector(), bottomLayer + size, nx, ny);
  };

  // top copy layer
  if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT)
  {
    dim3 dimBlock(TILE_SIZE, 1);
    dim3 dimGrid(nx / TILE_SIZE, 1);
    kernelTopCopyLayer<<<dimGrid, dimBlock>>>(
        h.elemVector(), hu.elemVector(), hv.elemVector(), topLayer + size, nx, ny);
  };

	cudaDeviceSynchronize();

}

void SWE_BlockCudaAwareMPI::synchGhostLayerAfterWrite()
{
#ifdef DBG
  cout << "Unpacking data ..." << flush << endl;
#endif

  // bottom copy layer
  if (boundary[BND_BOTTOM] == PASSIVE || boundary[BND_BOTTOM] == CONNECT)
  {
    dim3 dimBlock(TILE_SIZE, 1);
    dim3 dimGrid(nx / TILE_SIZE, 1);
    kernelBottomGhostBoundary<<<dimGrid, dimBlock>>>(
        h.elemVector(), hu.elemVector(), hv.elemVector(), bottomLayer, nx, ny);
  };

  // top copy layer
  if (boundary[BND_TOP] == PASSIVE || boundary[BND_TOP] == CONNECT)
  {
    dim3 dimBlock(TILE_SIZE, 1);
    dim3 dimGrid(nx / TILE_SIZE, 1);
    kernelTopGhostBoundary<<<dimGrid, dimBlock>>>(
        h.elemVector(), hu.elemVector(), hv.elemVector(), topLayer, nx, ny);
  };
	cudaDeviceSynchronize();

}
#endif
/**
 * register the row or column layer next to a boundary as a "copy layer",
 * from which values will be copied into the ghost layer or a neighbour;
 * @return	a SWE_Block1D object that contains row variables h, hu, and hv
 */
SWE_Block1D* SWE_BlockCudaAwareMPI::registerCopyLayer(BoundaryEdge edge){

  // for TOP and BOTTOM layer, the implementation is identical to that in SWE_Block
  // for LEFT and RIGHT layer, separate layers are used that avoid strided copies
  // when transferring memory between host and device memory
  switch (edge) {
    case BND_LEFT:
      return new SWE_Block1D( h.getColProxy(1), hu.getColProxy(1), hv.getColProxy(1));
    case BND_RIGHT:
      return new SWE_Block1D( h.getColProxy(nx), hu.getColProxy(nx), hv.getColProxy(nx));
#ifdef PACKING_DATA
    case BND_BOTTOM:
      synchCopyLayerBeforeRead();
      return bottomCopyLayer; //PHUONG: Here need to call synchronize so that it has data before first communication
    case BND_TOP:
      synchCopyLayerBeforeRead();
      return topCopyLayer;
#else
    case BND_BOTTOM:
      return new SWE_Block1D( h.getRowProxy(1), hu.getRowProxy(1), hv.getRowProxy(1));
    case BND_TOP:
      return new SWE_Block1D( h.getRowProxy(ny), hu.getRowProxy(ny), hv.getRowProxy(ny));
#endif
  };
  return NULL;
}

/**
 * "grab" the ghost layer at the specific boundary in order to set boundary values
 * in this ghost layer externally.
 * The boundary conditions at the respective ghost layer is set to PASSIVE,
 * such that the grabbing program component is responsible to provide correct
 * values in the ghost layer, for example by receiving data from a remote
 * copy layer via MPI communication.
 * @param	specified edge
 * @return	a SWE_Block1D object that contains row variables h, hu, and hv
 */
SWE_Block1D* SWE_BlockCudaAwareMPI::grabGhostLayer(BoundaryEdge edge){

  // for TOP and BOTTOM layer, the implementation is identical to that in SWE_Block
  // for LEFT and RIGHT layer, separate layers are used that avoid strided copies
  // when transferring memory between host and device memory
  boundary[edge] = PASSIVE;
  switch (edge) {
    case BND_LEFT:
      return new SWE_Block1D( h.getColProxy(0), hu.getColProxy(0), hv.getColProxy(0));
    case BND_RIGHT:
      return new SWE_Block1D( h.getColProxy(nx+1), hu.getColProxy(nx+1), hv.getColProxy(nx+1));
#ifdef PACKING_DATA
    case BND_BOTTOM:
      return bottomGhostLayer;
    case BND_TOP:
      return topGhostLayer;
#else
    case BND_BOTTOM:
      return new SWE_Block1D( h.getRowProxy(0), hu.getRowProxy(0), hv.getRowProxy(0));
    case BND_TOP:
      return new SWE_Block1D( h.getRowProxy(ny+1), hu.getRowProxy(ny+1), hv.getRowProxy(ny+1));
#endif
  };
  return NULL;
}

/**
 * Print some available information about the CUDA devices.
 */
void SWE_BlockCudaAwareMPI::printDeviceInformation()
{
	tools::Logger::logger.printString("Printing device information");

  //! id of the CUDA device.
  int l_deviceId;
  cudaGetDevice(&l_deviceId);

  //! total number of CUDA devices on this host.
  int l_deviceCount;
  cudaGetDeviceCount(&l_deviceCount);

  //! drive and runtime version
  int l_driverVersion, l_runtimeVersion;
  cudaDriverGetVersion(&l_driverVersion);
  cudaRuntimeGetVersion(&l_runtimeVersion);

  //! device properties
  cudaDeviceProp l_deviceProperty;
  cudaGetDeviceProperties(&l_deviceProperty, l_deviceId);

  // print information about the current device

  tools::Logger::logger.cout() << "Current CUDA device (relative to host): " << l_deviceId
                     << " ( " << l_deviceCount << " in total)" << std::endl;

  tools::Logger::logger.cout() << "CUDA device properties: "
                     << l_deviceProperty.name << " (name), "
                     << l_driverVersion << "/" << l_runtimeVersion << " (driver/runtime version), "
                     << l_deviceProperty.major << "." << l_deviceProperty.minor << " (compute capability)"
                     << std::endl;
}

void SWE_BlockCudaAwareMPI::init(int i_cudaDevice)
{
	  tools::Logger::logger.setProcessRank(i_cudaDevice);

	  cudaSetDevice(i_cudaDevice);

	  // check for a valid CUDA device id
	  #ifndef NDEBUG
	  int l_deviceCount;
	  cudaGetDeviceCount(&l_deviceCount);
	  assert( (i_cudaDevice >= 0) && (i_cudaDevice < l_deviceCount) );
	  #endif

	  printDeviceInformation();

	  // Make sure the cuda device is reset at exit
	  atexit( SWE_BlockCudaAwareMPI::finalize );
}

void SWE_BlockCudaAwareMPI::finalize()
{
	// reset the cuda device
	tools::Logger::logger.printString("Resetting the CUDA devices");
	cudaDeviceReset();
}
