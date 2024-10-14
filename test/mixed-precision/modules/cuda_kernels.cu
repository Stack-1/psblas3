#include <functional> // needed for bind
#include <stdio.h>

#define THREADS_PER_BLOCK 512


extern "C"  {void single_to_double(float *vect_single, double *vect_double, int size);}
extern "C"  {void d_to_r_kernel(double *r_double,float *r_single,double *d_double,float *d_single,int size);}
extern "C"  {void geaxpby_double_to_single(float *vect_single, double *vect_double, float *second_vect_single, int size, float alpha);}
extern "C"  {void geaxpby_double_to_double(float *vect_single, double *vect_double, int size, float beta);}
extern "C"  {void geaxpby_single(float *vect_single, float *second_vect_single, int size, float alpha);}

__global__ void convert_single_to_double(float *vect_single_dev, double *vect_double_dev, int size)
{
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid < size){
        vect_double_dev[tid] = static_cast<double>(vect_single_dev[tid]);
    }

}


void single_to_double(float *vect_single, double *vect_double, int size){
    float *vect_single_dev;
    double *vect_double_dev;
    int number_of_blocks;
    cudaError_t err;

    cudaMalloc((void**)&vect_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&vect_double_dev, sizeof(double)*size);

    cudaMemcpy(vect_single_dev, vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);

    number_of_blocks = size / THREADS_PER_BLOCK;
    if(number_of_blocks == 0){
        number_of_blocks = 1;
    }

    convert_single_to_double<<< number_of_blocks , THREADS_PER_BLOCK>>>(vect_single_dev,vect_double_dev,size);
    
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    
    if (err != cudaSuccess) 
        printf("Error in single_to_double: %s\n", cudaGetErrorString(err));


    cudaMemcpy(vect_double, vect_double_dev, sizeof(double) * size, cudaMemcpyDeviceToHost);
    
    cudaFree(vect_single_dev);
    cudaFree(vect_double_dev);
}




__global__ void d_to_r(double *r_double_dev,float *r_single_dev,double *d_double_dev,float *d_single_dev,int size)
{
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid < size){
        r_single_dev[tid]           = static_cast<float>(r_double_dev[tid]);
        d_double_dev[tid]           = r_double_dev[tid];
        d_single_dev[tid]           = static_cast<float>(d_double_dev[tid]);
    }

}


void d_to_r_kernel(double *r_double,float *r_single,double *d_double,float *d_single,int size){
    float *r_single_dev, *d_single_dev;
    double *r_double_dev, * d_double_dev;
    int number_of_blocks;
    cudaError_t err;

    cudaMalloc((void**)&r_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&d_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&r_double_dev, sizeof(double)*size);
    cudaMalloc((void**)&d_double_dev, sizeof(double)*size);

    cudaMemcpy(r_double_dev, r_double, sizeof(double) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_double_dev, d_double, sizeof(float) * size, cudaMemcpyHostToDevice);

    number_of_blocks = size / THREADS_PER_BLOCK;
    if(number_of_blocks == 0){
        number_of_blocks = 1;
    }

    d_to_r<<< number_of_blocks , THREADS_PER_BLOCK>>>(r_double_dev,r_single_dev,d_double_dev,d_single_dev,size);
    
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    
    if (err != cudaSuccess) 
        printf("Error in d_to_r_kernel: %s\n", cudaGetErrorString(err));

    cudaMemcpy(r_single, r_single_dev, sizeof(float) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(d_double, d_double_dev, sizeof(double) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(d_single, d_single_dev, sizeof(float) * size, cudaMemcpyDeviceToHost);


    
    cudaFree(r_single_dev);
    cudaFree(d_double_dev);
    cudaFree(d_single_dev);
    cudaFree(r_double_dev);
}

__global__ void geaxpby_s_d(float *vect_single_dev, double *vect_double_dev, float *second_vect_single_dev, float alpha, int size)
{
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid < size){
        vect_double_dev[tid]            = static_cast<double>(vect_single_dev[tid]) - ( static_cast<double>(alpha) * static_cast<double>(second_vect_single_dev[tid]) );
        vect_single_dev[tid]            = static_cast<float>(vect_double_dev[tid]);
    }

}


void geaxpby_double_to_single(float *vect_single, double *vect_double, float *second_vect_single, int size, float alpha){
    float *vect_single_dev, *second_vect_single_dev;
    double *vect_double_dev;
    int number_of_blocks;
    cudaError_t err;

    cudaMalloc((void**)&vect_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&second_vect_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&vect_double_dev, sizeof(double)*size);

    cudaMemcpy(vect_single_dev, vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(second_vect_single_dev, second_vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);

    number_of_blocks = size / THREADS_PER_BLOCK;
    if(number_of_blocks == 0){
        number_of_blocks = 1;
    }


    geaxpby_s_d<<< number_of_blocks , THREADS_PER_BLOCK>>>(vect_single_dev,vect_double_dev,second_vect_single_dev, alpha, size);
    
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    
    if (err != cudaSuccess) 
        printf("Error in geaxpby_double_to_single: %s\n", cudaGetErrorString(err));


    cudaMemcpy(vect_double, vect_double_dev, sizeof(double) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(vect_single, vect_single_dev, sizeof(float) * size, cudaMemcpyDeviceToHost);

    
    cudaFree(vect_single_dev);
    cudaFree(vect_double_dev);
    cudaFree(second_vect_single_dev);
}



__global__ void geaxpby_d_d(float *vect_single_dev, double *vect_double_dev, float *second_vect_single_dev, float beta, int size)
{
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid < size){
        vect_double_dev[tid]                    = static_cast<double>(vect_single_dev[tid]) + ( static_cast<double>(beta) * vect_double_dev[tid] );
        second_vect_single_dev[tid]             = static_cast<float>(vect_double_dev[tid]);
    }

}


void geaxpby_double_to_double(float *vect_single, double *vect_double, float *second_vect_single , int size, float beta){
    float *vect_single_dev, *second_vect_single_dev;
    double *vect_double_dev;
    int number_of_blocks;
    cudaError_t err;

    cudaMalloc((void**)&vect_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&second_vect_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&vect_double_dev, sizeof(double)*size);

    cudaMemcpy(vect_single_dev, vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(second_vect_single_dev, second_vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(vect_double_dev, vect_double, sizeof(double) * size, cudaMemcpyHostToDevice);

    number_of_blocks = size / THREADS_PER_BLOCK;
    if(number_of_blocks == 0){
        number_of_blocks = 1;
    }


    geaxpby_d_d<<< number_of_blocks , THREADS_PER_BLOCK>>>(vect_single_dev,vect_double_dev, second_vect_single_dev, beta, size);
    
    err = cudaGetLastError();
    
    if (err != cudaSuccess) 
        printf("Error in geaxpby_double_to_single: %s\n", cudaGetErrorString(err));


    cudaMemcpy(vect_double, vect_double_dev, sizeof(double) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(second_vect_single, second_vect_single_dev, sizeof(float) * size, cudaMemcpyDeviceToHost);

    
    cudaFree(vect_single_dev);
    cudaFree(vect_double_dev);
    cudaFree(second_vect_single_dev);
}




__global__ void geaxpby_s(float *vect_single_dev, float *second_vect_single_dev, float alpha, int size)
{
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid < size){
        second_vect_single_dev[tid] = second_vect_single_dev[tid] + ( alpha * vect_single_dev[tid] );
    }

}


void geaxpby_single(float *vect_single, float *second_vect_single, int size, float alpha){
    float *vect_single_dev, *second_vect_single_dev;
    int number_of_blocks;
    cudaError_t err;

    cudaMalloc((void**)&vect_single_dev, sizeof(float)*size);
    cudaMalloc((void**)&second_vect_single_dev, sizeof(float)*size);

    cudaMemcpy(vect_single_dev, vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(second_vect_single_dev, second_vect_single, sizeof(float) * size, cudaMemcpyHostToDevice);

    number_of_blocks = size / THREADS_PER_BLOCK;
    if(number_of_blocks == 0){
        number_of_blocks = 1;
    }


    geaxpby_s<<< number_of_blocks , THREADS_PER_BLOCK>>>(vect_single_dev,second_vect_single_dev, alpha, size);
    
    err = cudaGetLastError();
    
    if (err != cudaSuccess) 
        printf("Error in geaxpby_single: %s\n", cudaGetErrorString(err));


    cudaMemcpy(second_vect_single, second_vect_single_dev, sizeof(float) * size, cudaMemcpyDeviceToHost);

    
    cudaFree(vect_single_dev);
    cudaFree(second_vect_single_dev);
}
