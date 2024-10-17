/*
 * spGPU - Sparse matrices on GPU library.
 *  
 * Copyright (C) 2010 - 2015
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
#include "cudadebug.h"

#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#undef GEN_SPGPU_HELL_NAME
#undef X_TEX
#define X_TEX CONCAT(x_tex_, FUNC_SUFFIX)

__device__ __host__ static float zero_float() { return 0.0f; }
__device__ __host__ static bool float_isNotZero(float x) { return x != 0.0f; }

__device__ static double mixed_fma(float a, double b, float c) { return static_cast<double>(a)*b + static_cast<double>(c); }
__device__ static double mixed_add(float a, double b) { return static_cast<double>(a)+b; }
__device__ static double mixed_mul(float a, double b) { return static_cast<double>(a)*b; }

__device__ static float float_fma(float a, float b, float c) { return a*b + c; }
__device__ static float float_add(float a, float b) { return a+b; }
__device__ static float float_mul(float a, float b) { return a*b; }


__device__ static float readValue_float(float fetch) { return fetch; }

// host or c.c >= 1.3 
#if (__CUDA_ARCH__ >= 130) || (!__CUDA_ARCH__)
__device__ __host__ static double zero_double() { return 0.0; }
__device__ __host__ static bool double_isNotZero(double x) { return x != 0.0; }


__device__ static double readValue_double(int2 fetch) { return __hiloint2double (fetch.y, fetch.x); }
#endif
#if 0 
// Texture cache management
texture < TEX_FETCH_TYPE, 1, cudaReadModeElementType > X_TEX;

#define bind_tex_x(x) cudaBindTexture(NULL, X_TEX, x)
#define unbind_tex_x(x) cudaUnbindTexture(X_TEX)

__device__ static VALUE_TYPE 
fetchTex (int pointer)
{
	TEX_FETCH_TYPE fetch = tex1Dfetch (X_TEX, pointer);
	return CONCAT(readValue_,VALUE_TYPE) (fetch);
}
#endif
#if __CUDA_ARCH__ < 300
extern __shared__ int dynShrMem[]; 
#endif


#include "hell_spmv_base_template_mx.cuh"

void spgpuShellspmv_mx (spgpuHandle_t handle, 
	VALUE_TYPE_2* z, 
	const VALUE_TYPE_2 *y, 
	VALUE_TYPE alpha, 
	const VALUE_TYPE* cM, 
	const int* rP, 
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int avgNnzPerRow,	
	int rows, 
	const VALUE_TYPE *x, 
	VALUE_TYPE beta, 
	int baseIndex)
{

	int maxNForACall = max(handle->maxGridSizeX, THREAD_BLOCK*handle->maxGridSizeX);

	// maxNForACall should be a multiple of hackSize
	maxNForACall = (maxNForACall/hackSize)*hackSize;
	
	while (rows > maxNForACall) //managing large vectors
	{
	  
	  spgpuShellspmv_vanilla_mx (handle, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rIdx, avgNnzPerRow, maxNForACall, x, beta, baseIndex);

	  y = y + maxNForACall;
	  z = z + maxNForACall;
	  hackOffsets = hackOffsets + maxNForACall/hackSize;
	  rS = rS + maxNForACall;
	  
	  rows -= maxNForACall;
	}
	spgpuShellspmv_vanilla_mx (handle, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rIdx, avgNnzPerRow, rows, x, beta, baseIndex);

	cudaCheckError("CUDA error on hell_spmv");
}

