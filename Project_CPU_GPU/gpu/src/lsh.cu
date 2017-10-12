/**
 *	Locality senstive Hashing CPU code
 *  Implementation of routine function of LSH
 *
 *	@created 
 *	@Arvind and Divyanshu
 *
 *	The MIT License (MIT)
 *
 *	Copyright (c) 2017 Arvind and Divyanshu
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
 *	associated documentation files (the "Software"), to deal in the Software without restriction, 
 *	including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 *	and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 *	subject to the following conditions:
 *
 *		The above copyright notice and this permission notice shall be included in all copies or 
 *		substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
 *	LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 *	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
 *	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
 *	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <array>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "lsh.h"
#include <curand.h>
#include <curand_kernel.h>
#include <sys/time.h>
#include <iostream>
#include <cstring>

using namespace std;

#define BLOCK_SIZE 16

float **data;
table a[4096];
long get_usecss (void)
{
   struct timeval t;
   gettimeofday(&t,NULL);
   return t.tv_sec*1000000+t.tv_usec;
}
//utility functions
__device__ char * my_strcpy(char *dest, const char src){
  int i = 0;
  do {
    dest[i] = src;
  }
  while (src!= 0);
  return dest;
}
//utility functions
__device__ char * my_strcat(char *dest, const char src){
  int i = 0;
  while (dest[i] != 0) i++;
  my_strcpy(dest+i, src);
  return dest;
}


// Kernel Function for hyperplane
__global__ void cuda_hyperplane_kernel(float *d_data)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int seed = time_t(NULL);
	curandState s;
	curand_init(seed, blockIdx.x, 0, &s);
	d_data[i] = curand_normal(&s);
}
// Kernel Function for hash_matrix without shared memory concept
__global__ void cuda_hashmatrix_kernel(float *A,float *B,float *C,int d_p, int d_col,int d_rows)
{
	float result = 0;
    int Row = blockIdx.y*BLOCK_SIZE + threadIdx.y;
    int Col = blockIdx.x*BLOCK_SIZE + threadIdx.x;
    
    __shared__ float As[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ float Bs[BLOCK_SIZE][BLOCK_SIZE];

    for (int k = 0; k < (BLOCK_SIZE + d_rows - 1)/BLOCK_SIZE; k++) {
    	if (k*BLOCK_SIZE + threadIdx.x < d_rows && Row < d_p)
             As[threadIdx.y][threadIdx.x] = A[Row*d_rows + k*BLOCK_SIZE + threadIdx.x];
        else
             As[threadIdx.y][threadIdx.x] = 0.0;

        if (k*BLOCK_SIZE + threadIdx.y < d_rows && Col < d_col)
             Bs[threadIdx.y][threadIdx.x] = B[(k*BLOCK_SIZE + threadIdx.y)*d_col + Col];
        else
             Bs[threadIdx.y][threadIdx.x] = 0.0;

         __syncthreads(); 

        for (int n = 0; n < BLOCK_SIZE; ++n) 
            // if ((k*BLOCK_SIZE + n < d_rows && Row < d_p) && (k*BLOCK_SIZE + n < d_rows && Col < d_col))
                // result += A[Row*d_rows + k*BLOCK_SIZE + n] * B[(k*BLOCK_SIZE + n)*d_col + Col];
        	result+=As[threadIdx.y][n] * Bs[n][threadIdx.x];

    }

    if (Row < d_p && Col < d_col)
    {
    	if (result >=0){
			C[((blockIdx.y * blockDim.y + threadIdx.y)*d_col)+(blockIdx.x*blockDim.x)+threadIdx.x]=1;		
    	}
    	else
    	{
			C[((blockIdx.y * blockDim.y + threadIdx.y)*d_col)+(blockIdx.x*blockDim.x)+threadIdx.x]=0;
    	}
    
    }     	
}
__global__ void cuda_hashmatrix_merge_kernel(char **A ,table * B,int d_p,int d_col,int buckets)
{
	int p = 25;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
    // B[0].key = A[0];
    memcpy(B[0].key, A[0], sizeof(char)*p);
    B[0].count += 1;
    B[0].values[B[0].count]=0;
    __syncthreads();
    if(B[i].key == A[i])
    {
    	B[i].values[B[i].count]=i;
    	B[i].count += 1;
    }
    else
    {
    	memcpy(B[i].key, A[0], sizeof(char)*p);
        B[i].values[B[i].count] = i;
        B[i].count += 1;
    }
}
__global__ void cuda_cosine_Kernel (float *d_dataset,float *d_query,int *d_samples ,float *d_cosine_values, int rows, int col) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int index = d_samples[i];
	float dot = 0.0, denom_a = 0.0, denom_b = 0.0 ;
	for(unsigned int j = 0u; j < rows; j++) {
        dot += d_dataset[col * j + index] * d_query[1 * j ] ;
        denom_a += d_dataset[col * j + index] * d_dataset[col * j + index] ;
        denom_b += d_query[1 * j ] * d_query[1 * j ] ;
    }
    d_cosine_values[i] = dot / (sqrt(denom_a) * sqrt(denom_b));
    // dot = 0.0, denom_a = 0.0, denom_b = 0.0 ;
    //printf("%f\n", d_cosine_values[i]);
}

__global__ void cuda_hamming_Kernel(table * device_table, char * device_query, int * device_hamming_dist, int p) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int count = 0;
    
    for (int j = 0; j<p; j++) {
        if (device_table[i].key[j] != device_query[j]) {
            count++;
        }
    }
    device_hamming_dist[i] = count;
}

/*	Read input data set and return it as matrix*/
float ** lsh::read_data(std::string filename)
{
	string line;
	int row=0;
	int col=0;	
	data = new float*[58307];
	for(int i = 0; i < 58307; ++i)
    	data[i] = new float[4096];
	
	ifstream myfile (filename);
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			std::istringstream iss(line);
			int val;
			col=0;
			while ( iss >> val)
			{
			
				data[row][col]=val;
				col++;
			}
	  		row++;
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	// cout << row << endl << col;
	return data;
}

/*	Hyperplane matrix generation algorithm and return the hyperplane matrix [TODO number generated by random algorithm is not in range of 0 to 1] */
float * lsh::hyper_plane(int p, int rows)
{
	float *data;
	data = new float[p*rows];
    // define device variables
    float* d_data;
    // allocate the memory for the deivce varaible
    cudaMalloc((void **)&d_data,p*rows*sizeof(float));
    // calling the kernel function
    long start = get_usecss();
    cuda_hyperplane_kernel<<<p*rows/1024,1024>>>(d_data);
	long end = get_usecss();
	double dur = ((double)(end-start))/1000000;
	printf("Hyperplance kernel Compute Time = %f\n",dur);
    // copy back the result and send the output to cpu
	cudaMemcpy(data,d_data,p*rows*sizeof(float),cudaMemcpyDeviceToHost);	
	return data;
}

float * lsh::hash_matrix(float *dataset , float *hyperplane,int p ,int col, int rows)
{
	cudaError_t cudaStatus;
	// host varaibles
	float *data;
	data = new float[p*col];
	// device varaibles 
	float *d_dataset,*d_hyperplane,* d_data = 0;
	// allocate storage for the device
	cudaMalloc((void**)&d_hyperplane, sizeof(float) * p * rows);
	cudaMalloc((void**)&d_dataset, sizeof(float) * rows * col);
	cudaMalloc((void**)&d_data, sizeof(float)* p * col);	
	// cudaMemset(d_data, 0, sizeof(float)* p * col);
	// copy input to the device
	cudaStatus=cudaMemcpy(d_dataset, dataset, sizeof(float) * rows * col, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed! %s",cudaGetErrorString(cudaGetLastError()));
	}
	cudaStatus=cudaMemcpy(d_hyperplane, hyperplane, sizeof(float) * p * rows, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed! %s",cudaGetErrorString(cudaGetLastError()));
	}
	// calling the kernel function
    long start = get_usecss();
    const dim3 block_size(BLOCK_SIZE,BLOCK_SIZE);
    const dim3 num_blocks((col + block_size.x - 1)/block_size.x, (p + block_size.y - 1)/block_size.y);
    cuda_hashmatrix_kernel<<<num_blocks,block_size>>>(d_hyperplane,d_dataset,d_data, p, col,rows);
    cudaStatus = cudaMemcpy(data, d_data ,sizeof(float)*p*col,cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed! %s",cudaGetErrorString(cudaGetLastError()));
	}	
	long end = get_usecss();
	double dur = ((double)(end-start))/1000000;
	printf("Hash matrix kernel Compute Time = %f\n",dur);
    // copy back the result and send the output to cpu
	return data;
}

table * lsh::hashtable(float *dataset,int p ,int col, int rows , int *no_of_buckets)
{
	// cudaError_t cudaStatus;
	// cudaError_t err = cudaSuccess;
	// host varaibles
	// char *data[col];
	// char (*data)[col] = new char[col][25];
	// char ** data = new char*[col];
	char sum[p];
	// for (int i=0;i<col;i++) data[i] = new char[25];
	// data = new char[col];
	
	int buckets=1;
	int flag=0;	
	for (int i = 0; i < col; i++)
	{	
		for (int j = 0; j < p; j++)
		{
			if (dataset[col * j + i]>0) 
				sum[j] = '1';
			else
				sum[j] = '0';
		}
		
		if(i==0)
		{
			strcpy(a[i].key,sum);
			a[i].values[a[i].count]=i;
			a[i].count++;
			// buckets++;
		}
		else
		{
			flag=0;
			int k;
			for(k=0;k<buckets;k++)
			{
				if(strcmp(a[k].key,sum)==0)
				{
					a[k].values[a[k].count]=i;
					a[k].count++;
					flag=1;
				}
			}
			if(flag==0)
			{
				strcpy(a[k].key,sum);
				a[k].values[a[k].count]=i;
				a[k].count++;
				buckets++;
			}
		}
	}
	int total_counts = 0;
 	for (int i = 0; i< buckets; i++) {
 		// printf("COUNT = %d\n", a[i].count);
 		total_counts=total_counts+a[i].count;
 	}
 	*no_of_buckets=buckets;
	// // device varaibles 
	// char **d_dataset;
	// table *d_data;
	// // allocate storage for the device
	// cudaMalloc((void**)&d_dataset, sizeof(char)*col*25);
	// err = cudaMalloc((void**)&d_data, sizeof(table)*col);
 //    if (err != cudaSuccess)
 //    {
 //        fprintf(stderr, "Failed to allocate hash table (error code %s)!\n", cudaGetErrorString(err));
 //        exit(EXIT_FAILURE);
 //    }	
	// // copy input to the device
	// cudaStatus=cudaMemcpy(d_dataset, data, sizeof(char) * col, cudaMemcpyHostToDevice);
	// if (cudaStatus != cudaSuccess) {
	// 	fprintf(stderr, "cudaMemcpy1 failed! %s",cudaGetErrorString(cudaGetLastError()));
	// }
	// // calling the kernel function TODO implement correctly the kernel function
 //    long start = get_usecss();
 //    // const dim3 block_size(BLOCK_SIZE,BLOCK_SIZE);
 //    // const dim3 num_blocks((col + block_size.x - 1)/block_size.x, (p + block_size.y - 1)/block_size.y);
 //    cuda_hashmatrix_merge_kernel<<<col/1024,1024>>>(d_dataset,d_data, p, col,buckets);
 //    cudaStatus = cudaMemcpy(a, d_data,sizeof(table)*col,cudaMemcpyDeviceToHost);
	// if (cudaStatus != cudaSuccess) {
	// 	fprintf(stderr, "cudaMemcpy2 failed! %s",cudaGetErrorString(cudaGetLastError()));
	// }
	// cout << a[0].key <<"\t" << a[1].key<<"\t" << a[2].key<<"\t"<< a[3].key<<endl;
	// long end = get_usecss();
	// double dur = ((double)(end-start))/1000000;
	// printf("Hash matrix merge kernel Compute Time = %f\n",dur);
 //    // copy back the result and send the output to cpu
	return a;	
}
int * lsh::hamming_distance(table *hash , char query[],int p, int buckets,int col)
{
	// cudaError_t err = cudaSuccess;
    int * host_hamming_dist;
    host_hamming_dist = (int *) malloc(buckets*sizeof(int));
    //  device variables
    char * device_query;
    table * device_table;
    int * device_hamming_dist;
    // allocating the memory
    cudaMalloc((void **) &device_table, col*sizeof(table));
    cudaMalloc((void **) &device_query, p*sizeof(char));
    cudaMalloc((void **) &device_hamming_dist, buckets*sizeof(int));
    // Copy data to device memory
    cudaMemcpy(device_table, hash, col*sizeof(table), cudaMemcpyHostToDevice);
    cudaMemcpy(device_query, query, p*sizeof(char), cudaMemcpyHostToDevice);
    // calling the kernel
    cuda_hamming_Kernel<<<1, buckets>>> (device_table, device_query, device_hamming_dist, p);
    // copy result back
    cudaMemcpy(host_hamming_dist, device_hamming_dist, buckets*sizeof(int), cudaMemcpyDeviceToHost);

    return host_hamming_dist;
}
float * lsh::cosine_distance(float *dataset,float *query ,int *samples, int rows,int col,int count)
{
	// cudaError_t err = cudaSuccess;
    float * host_cosine_dist;
    host_cosine_dist = (float *) malloc(col*sizeof(float));
    //  device variables
    float * device_query;
    float * device_data;
    float * device_cosine_dist;
	int *device_samples;
	// printf("%f\n",dataset[0]);
	// printf("%d\n",count);
	// printf("%f\n",query[2]);
	// printf("%d\n",samples[0]);
    // allocating the memory
    cudaMalloc((void **) &device_data, col*rows*sizeof(float));
    cudaMalloc((void **) &device_query, rows*sizeof(float));
    cudaMalloc((void **) &device_samples, col*sizeof(int));
    cudaMalloc((void **) &device_cosine_dist, col*sizeof(float));
    // Copy data to device memory
    cudaMemcpy(device_data,dataset, col*rows*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(device_query, query, rows*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(device_samples, samples, col*sizeof(int), cudaMemcpyHostToDevice);
    // calling the kernel
    if (count>1024)
	    cuda_cosine_Kernel<<<1, 1024>>> (device_data,device_query,device_samples , device_cosine_dist, rows,col);
    else	
    	cuda_cosine_Kernel<<<1, count>>> (device_data,device_query,device_samples ,device_cosine_dist, rows,col);
    // copy result back
    cudaMemcpy(host_cosine_dist, device_cosine_dist, col*sizeof(float), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < 10; ++i)
    // {
    // 	 code 
    // 	printf("%f\n", host_cosine_dist[i]);
    // }
    return host_cosine_dist;
}