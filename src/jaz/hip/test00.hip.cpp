#include <hip/hip_runtime.h>
#include "test00.h"
#include "kernels/add.h"
#include <iostream>

void HipTest00 :: run()
{
  int N = 1<<20; // 1M elements

  // Allocate Unified Memory -- accessible from CPU or GPU
  float *x, *y;
  hipMallocManaged(&x, N*sizeof(float));
  hipMallocManaged(&y, N*sizeof(float));

  // initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
	x[i] = 1.0f;
	y[i] = 2.0f;
  }

  // Run kernel on 1M elements on the GPU
  hipLaunchKernelGGL(add, 1, 1, 0, 0, N, x, y);

  // Wait for GPU to finish before accessing on host
  hipDeviceSynchronize();

  // Check for errors (all values should be 3.0f)
  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
  {
	maxError = fmax(maxError, fabs(y[i]-3.0f));
  }

  std::cout << "Max error: " << maxError << std::endl;

  // Free memory
  hipFree(x);
  hipFree(y);
}
