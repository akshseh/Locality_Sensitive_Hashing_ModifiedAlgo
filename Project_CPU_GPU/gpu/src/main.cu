/**
 *	Locality senstive Hashing CPU code
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
// #include <sstream>
// #include <random>
#include <iostream>
#include <string.h>
#include <array>
#include <unordered_map>
#include <string>
#include <vector>
#include "lsh.h"
#include <sys/time.h>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <map>
#include <thrust/version.h>
#include <thrust/sort.h>
using namespace std;

long get_usecs (void)
{
   struct timeval t;
   gettimeofday(&t,NULL);
   return t.tv_sec*1000000+t.tv_usec;
}

double cosine_similarity(float *A, float *B, unsigned int Vector_Length)
{
    double dot = 0.0, denom_a = 0.0, denom_b = 0.0 ;
 	for(unsigned int i = 0u; i < Vector_Length; ++i) {
        dot += A[i] * B[i] ;
        denom_a += A[i] * A[i] ;
        denom_b += B[i] * B[i] ;
    }
    return dot / (sqrt(denom_a) * sqrt(denom_b)) ;
}

int main()
{
	// string path = "../../data/test";
	int p=25,rows=58307,col=512;
	cout << "Enter the value of P \t";
	cin >> p;
	cout << "Enter the value of col as per dataset \t";
	cin >> col;
	int knn=100;
	cout << "Enter the value of k \t";
	cin >> knn;

	long start = get_usecs();
	float * hyperplane= lsh::hyper_plane(p,rows);
	long end = get_usecs();
	double dur = ((double)(end-start))/1000000;
	printf("Hyperplane with memory transfer Time = %f\n",dur);

	
	long tstart = get_usecs();
	string path = "../data/dataset_1024.txt";
	start = get_usecs();
	
	


	float ** dataset=lsh::read_data(path);
	end = get_usecs();
	dur = ((double)(end-start))/1000000;
	float  *datasets = new float[rows*col];
	float ** trans_data;
	trans_data = new float*[col];
	for(int i = 0; i < col; ++i)
    	trans_data[i] = new float[58307];
	
	for (int h = 0; h < rows; h++){
	    for (int w = 0; w < col; w++){
	        datasets[col * h + w] = dataset[h][w];
	        trans_data[w][h]=dataset[h][w];
	    }

	}
	printf("Dataset Load Time = %f\n",dur);
	//  TODO no rows are hardcorded 
	// float * trans_data = new float[col*rows];
	// for(int i = 0; i < rows; ++i)
	// {
	// 	for(int j = 0; j < col; ++j)
 //        {
 //        	trans_data[rows * j + i]=dataset[i][j];
 //     	}
	// }

	
	cudaDeviceSynchronize();
 	start = get_usecs();
	float * hash_matrix = lsh::hash_matrix(datasets,hyperplane,p,col,rows);
	end = get_usecs();
	dur = ((double)(end-start))/1000000;
    printf("Hash matrix  Time = %f\n",dur);
    int no_of_buckets=0;
    start = get_usecs();
    table * hasmerge = lsh::hashtable(hash_matrix,p ,col,rows,&no_of_buckets);
    end = get_usecs();
	dur = ((double)(end-start))/1000000;
	printf("Hash table  Time = %f\n",dur);
	start = get_usecs();
	string path1 = "../data/query1.txt";
    float ** query=lsh::read_data(path1);
    float query_trans[1][rows];
	float  *tquery= new float[rows];
	float  *nquery= new float[rows];
    for (int i = 0; i < rows; i++)
    {
    	for (int j = 0; j < 1; j++)
    	{
    		query_trans[j][i]=query[i][j];
	        tquery[rows * j + i] = query[i][j];
	        nquery[1* i + j] = query[i][j];
    	}
     	 
    }
	char host_query[p];
	for(int i = 0; i <p; i++)
	{
		for(int j = 0; j<1; j++)
		{
			float sum=0;
            for(int k = 0; k <rows; k++)
            {
                sum += hyperplane[col*i+k] * nquery[1*k+j]; // 
                if (sum >=0 ){
                	host_query[i] = '1';
                } 
                else{
                	host_query[i] = '0';
                }
            }
        }
	}
	// std::cout<<host_query<<endl;
	// Doing a K nearest search
	
	int *rank;
	int samples[col];
	float *cosine_rank;
	int cosine_values[col]; 
	int k=0;
	// searching for query key with the hashtable key
	for (int i = 0; i < no_of_buckets; ++i)
	{
		if(strcmp(hasmerge[i].key,host_query)!=0)
		{
			// no match of hash code
			rank=lsh::hamming_distance(hasmerge,host_query,p,no_of_buckets,col);
			
			for(i=0;i<no_of_buckets;i++)
			{
				if(rank[i] < 15)
				{
					for(int j=0;j<hasmerge[i].count;j++)
					{
						// cout <<hasmerge[i].values[j] <<"\t";
						samples[k]=hasmerge[i].values[j];
						k++;
					}
					// cout<<endl;
				}
			}
			// int count=0;
			// cosine_rank=lsh::cosine_distance(hasmerge,host_query,p,no_of_buckets,col);
			for (int j = 0; j < k; j++)
			{
				cosine_values[j]=samples[j];
				// cosine_rank[j]=cosine_similarity(trans_data[samples[j]],query_trans[0],rows);
				// ranks[]=samples[k];
				//cout << samples[k]<<"\t"<<cosine_similarity(trans_data[k],query_trans[0],rows)<<"\n";
			}
			cosine_rank=lsh::cosine_distance(datasets,tquery,samples,rows,col,k);
		    // for (int i = 0; i < 10; ++i)
		    // {
		    // 	 code 
		    // 	printf("%f\n", cosine_rank[i]);
		    // }
			thrust::sort_by_key(cosine_rank, cosine_rank + k, cosine_values,thrust::greater<float>());
			for(int i=0;i<knn;i++)
			{
	   			printf("Simalrity %f Sample column %i\n",cosine_rank[i],cosine_values[i]);	
			}
		}
		else
		{
			// int k=0;
			int temp = knn;
			for(int j=0;j<hasmerge[i].count;j++)
			{
				samples[k]=hasmerge[i].values[j];
				k++;
			}
	  		if(k==temp)
	  		{
	  			for (int j = 0; j < k; j++)
				{
					cosine_values[j]=samples[j];
					cosine_rank[j]=cosine_similarity(trans_data[samples[j]],query_trans[0],rows);
					// ranks[cosine_similarity(trans_data[k],query_trans[0],rows)]=samples[k];
				}

	  		}
	  		else
	  		{
	  			int match=i;
		  		rank=lsh::hamming_distance(hasmerge,host_query,p,no_of_buckets,col);
				// int k=0;
				for(i=0;i<no_of_buckets;i++)
				{
					if(i==match)
						continue;
					else if(rank[i] < 15)
					{
						for(int j=0;j<hasmerge[i].count;j++)
						{
							// cout <<hasmerge[i].values[j] <<"\t";
							samples[k]=hasmerge[i].values[j];
							k++;
						}
						// cout<<endl;
					}
				}
				// cosine_rank=lsh::cosine_distance(hasmerge,host_query,p,no_of_buckets,col);
				for (int j = 0; j < k; j++)
				{
					cosine_values[j]=samples[j];
					cosine_rank[j]=cosine_similarity(trans_data[samples[j]],query_trans[0],rows);
					// ranks[]=samples[k];
					//cout << samples[k]<<"\t"<<cosine_similarity(trans_data[k],query_trans[0],rows)<<"\n";
				}				
	  		}
	  		thrust::sort_by_key(cosine_rank, cosine_rank + k, cosine_values,thrust::greater<float>());
  			for(int i=0;i<knn;i++)
			{
	   			printf("Simalrity %f Sample column %i\n",cosine_rank[i],cosine_values[i]);	
			}
		}
	}
	end = get_usecs();
	dur = ((double)(end-start))/1000000;
	printf("\nQuery Search Time = %f\n",dur);
    long tend = get_usecs();
	double tdur = ((double)(tend-tstart))/1000000;
	printf("Total Program Time = %f\n",tdur);
  	return 0;
}
