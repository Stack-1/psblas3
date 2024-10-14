#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_atomic_functions.h"
#include <stdio.h>
#include <stdlib.h>
#define N 10
#define THREADS_PER_BLOCK 1

__global__ void dot(int *a, int *b, int *c)
{
    __shared__ int temp[THREADS_PER_BLOCK];
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    temp[threadIdx.x] = a[index] * b[index];

    __syncthreads();

    if (threadIdx.x == 0)
    {
        int sum = 0;
        for (int i = 0; i < THREADS_PER_BLOCK; i++)
        {
            sum += temp[i];
        }
        atomicAdd(c, sum);
    }
}

int main()
{
    int *a, *b, *c;
    int *dev_a, *dev_b, *dev_c;
    int size = N * sizeof(int);

    //allocate space for the variables on the device
    cudaMalloc((void **)&dev_a, size);
    cudaMalloc((void **)&dev_b, size);
    cudaMalloc((void **)&dev_c, sizeof(int));

    //allocate space for the variables on the host
    a = (int *)malloc(size);
    b = (int *)malloc(size);
    c = (int *)malloc(sizeof(int));

    //this is our ground truth
    int sumTest = 0;
    //generate numbers
    for (int i = 0; i < N; i++)
    {
        a[i] = rand() % 10;
        b[i] = rand() % 10;
        sumTest += a[i] * b[i];
    }
    printf("A: ");
    for(int i = 0; i < N; i++){
        printf("%d ",a[i]);
    }
    puts("");
    printf("B: ");
    for(int i = 0; i < N; i++){
        printf("%d ",b[i]);
    }
    puts("");
    *c = 0;

    cudaMemcpy(dev_a, a, sizeof(int) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, sizeof(int) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_c, c, sizeof(int), cudaMemcpyHostToDevice);

    dot<<< N / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> >(dev_a, dev_b, dev_c);
    cudaError_t err = cudaGetLastError();
    
    if (err != cudaSuccess) 
        printf("Error: %s\n", cudaGetErrorString(err));


    cudaMemcpy(c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);

    printf("GPU product: %d\n", *c);
    printf("CPU product: %d\n", sumTest);

}