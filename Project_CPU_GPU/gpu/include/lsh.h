/**
 *	lsh_cpu.H
 *
 */
#ifndef _lsh_H_INCLUDED
#define _lsh_H_INCLUDED
#include <cuda_runtime.h>
#include <string.h>
#include <iostream>
#include <unordered_map>
#include <array>
#include <vector>
#include <string>


typedef struct 
{
    char key[250];
    int  values[10710];
    int count = 0;
}table;

namespace lsh{
	/**
	 *	@param text - File path to read the file and return the matrix data
	 */
	float ** read_data(std::string path);
	/**
	 *	@param text - File path to read the file and return the matrix data
	 */
	float * hyper_plane(int p, int rows);	
	 
	float * hash_matrix(float *data , float *hyperplane, int p ,int n, int rows);

	table * hashtable(float *dataset,int p ,int col, int rows, int *no_buckets);

	int * hamming_distance(table *hash , char query[],int p , int buckets,int col);

	float * cosine_distance(float *dataset,float *query ,int *samples, int rows,int col,int count);
}

#endif
