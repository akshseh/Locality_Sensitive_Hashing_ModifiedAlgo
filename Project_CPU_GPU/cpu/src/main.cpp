// #include <sstream>
// #include <random>

// #include <iostream>
// #include <string.h>
// #include <array>
// #include <unordered_map>
// #include <string>
// #include <vector>
// #include "lsh.h"
// #include <sys/time.h>
// #include <math.h>
// #include <iterator>
// #include <algorithm>
// #include <map>
// using namespace std;
// // to do add it to the lsh header file Hamming distance help taken from stackoverflow.
// int GetHammingDistance(const std::string a, const std::string b)
// {
//     std::string::const_iterator a_it = a.begin();
//     std::string::const_iterator b_it = b.begin();
//     std::string::const_iterator a_end = a.end();
//     std::string::const_iterator b_end = b.end();
//     int value = 0;
//     while (a_it != a_end && b_it != b_end)
//     {
//         if (*a_it != *b_it){ 
//         	++value;
//         }
//         ++a_it;
//         ++b_it;
//     }
//     return value;
// }

// long get_usecs (void)
// {
//    struct timeval t;
//    gettimeofday(&t,NULL);
//    return t.tv_sec*1000000+t.tv_usec;
// }

// double cosine_similarity(float *A, float *B, unsigned int Vector_Length)
// {
//     double dot = 0.0, denom_a = 0.0, denom_b = 0.0 ;
//  	for(unsigned int i = 0u; i < Vector_Length; ++i) {
//         dot += A[i] * B[i] ;
//         denom_a += A[i] * A[i] ;
//         denom_b += B[i] * B[i] ;
//     }
//     return dot / (sqrt(denom_a) * sqrt(denom_b)) ;
// }

// int main()
// {
// 	// lsh::example_printf(10);
// 	// string path = "../../data/test";
// 	long tstart = get_usecs();
// 	string path = "../../data/dataset_2048.txt";
// 	int p=25,rows=58307,col=2048;
// 	int k = 40;
// 	long start = get_usecs();
// 	float ** dataset=lsh::read_data(path);
// 	long end = get_usecs();
// 	double dur = ((double)(end-start))/1000000;
// 	printf("Dataset Load Time = %f\n",dur);
// 	//  TODO no rows are hardcorded 
// 	float ** trans_data;
// 	trans_data = new float*[col];
// 	for(int i = 0; i < col; ++i)
//     	trans_data[i] = new float[58307];
	
// 	for(int i = 0; i < rows; ++i)
// 	{
// 		for(int j = 0; j < col; ++j)
//         {
//             trans_data[j][i]=dataset[i][j];
//      	}
// 	}
// 	 start = get_usecs();
// 	float ** hyperplane= lsh::hyper_plane(p,rows);    
// 	 end = get_usecs();
// 	 dur = ((double)(end-start))/1000000;
// 	printf("Hyperplane  Time = %f\n",dur);
// 	 start = get_usecs();
// 	float ** hash_matrix = lsh::hash_matrix(dataset,hyperplane,p,col,rows);
// 	 end = get_usecs();
// 	 dur = ((double)(end-start))/1000000;
//     printf("Hash matrix  Time = %f\n",dur);
// 	start = get_usecs();
// 	std::unordered_map<std::string, std::vector<int> > hash_tables = lsh::hash_table(hash_matrix,col,p); 
// 	end = get_usecs();
// 	dur = ((double)(end-start))/1000000;
// 	printf("Hash table  Time = %f\n",dur);
// 	std::unordered_map<std::string, std::vector<int> >:: iterator itr;
// 	start = get_usecs();
// 	string path1 = "../../data/query_738.txt";
//     float query_trans[1][rows];
//     float ** query=lsh::read_data(path1);
//     for (int i = 0; i < rows; i++)
//     {
//     	for (int j = 0; j < 1; j++)
//     	{
//     		query_trans[j][i]=query[i][j];
//     	}
     	 
//     }
// 	// generating the hash code for the query
//  	std::string code="";
// 	hash_matrix = lsh::hash_matrix(query,hyperplane,p,1,rows);
//   	for (int i = 0; i <p; i++)
// 	{
// 		//std::string code = "";
// 		for(int j=0; j< 1 ; j++)
// 		{
// 			code.append(std::to_string((int)hash_matrix[i][j]));
// 		}
// 	}
// 	// Doing a K nearest search based on string compare function TODO cosine etc
// 	int rank[hash_tables.size()];
// 	int samples[hash_tables.size()];
// 	std::map<float, int> ranks;
// 	std::string code1=code;
// 	// std::cout<<code1<<endl;
// 	// for (int i=0;i<rows;i++)
// 	// 	std::cout<<trans_data[0][i]<<"\t";
// 	// std::cout<<endl;
// 	// for (int i=0;i<rows;i++)
// 	// 	std::cout<<query_trans[0][i]<<"\t";
// 	// cout<<endl;
// 	// std::cout<<cosine_similarity(trans_data[0],query_trans[0],rows)<<endl;
// 	int j=0;
// 	if (hash_tables.find(code) == hash_tables.end()){
// 		// cout << "No exact match";
// 		int i=0;
// 		for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
// 		{
// 				rank[i]=GetHammingDistance(code1,itr->first);
// 				i++;
// 		}
// 		i=0;
// 		for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
// 		{
// 			if(rank[i] < 10)
// 			{
// 				// cout <<itr->first;
// 				std::vector<int> val = itr->second;
// 			    while (!val.empty()) {
// 		   			 auto it = val.rbegin();
// 		    		 // printf("%d\n", *it);
// 		    		 samples[j]=*it;
// 		    		 val.pop_back();
// 		    		 j++;
// 		  		}
// 			}
// 			i++;
// 		}
// 		for (int k = 0; k < j; k++)
// 		{
// 			ranks[cosine_similarity(trans_data[samples[k]],query_trans[0],rows)]=samples[k];
// 			//cout << samples[k]<<"\t"<<cosine_similarity(trans_data[k],query_trans[0],rows)<<"\n";
// 		}
//   		std::cout << "Nearest Samples to given query are:" <<endl;
//   		for (auto it =ranks.rbegin(); it !=ranks.rend(); it++) { // calls a_map.begin() and a_map.end()
// 			if(k>0){
// 				std::cout << k <<" Similarity "<<it->first <<" Sample column " <<it->second << '\n';	
// 				k--;
// 			}
    		
//   		}
// 	}
// 	// Key found part
// 	else
// 	{
// 		int temp = k;
// 		std::vector<int> val = hash_tables[code1];
// 	    while (!val.empty()) {
//    			 auto it = val.rbegin();
//     		 // printf("%d\n", *it);
//     		 samples[j]=*it;
//     		 val.pop_back();
//     		 j++;
//   		}
//   		if(j==temp)
//   		{
//   			for (int k = 0; k < j; k++)
// 			{
// 				ranks[cosine_similarity(trans_data[samples[k]],query_trans[0],rows)]=samples[k];
// 			}

//   		}
//   		else
//   		{
// 			int i=0;
// 			for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
// 			{
// 				if(itr->first == code1)
// 					continue;
// 				else{
// 					rank[i]=GetHammingDistance(code1,itr->first);
// 					i++;
// 				}
// 			}
// 			i=0;
// 			for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
// 			{
// 				if(rank[i] < 15)
// 				{
// 					// cout <<itr->first;
// 					std::vector<int> val = itr->second;
// 				    while (!val.empty()) {
// 			   			 auto it = val.rbegin();
// 			    		 // printf("%d\n", *it);
// 			    		 samples[j]=*it;
// 			    		 val.pop_back();
// 			    		 j++;
// 			  		}
// 				}
// 				i++;
// 			}
// 			for (int k = 0; k < j; k++)
// 			{
// 				ranks[cosine_similarity(trans_data[samples[k]],query_trans[0],rows)]=samples[k];
// 				//cout << samples[k]<<"\t"<<cosine_similarity(trans_data[k],query_trans[0],rows)<<"\n";
// 			}

//   		}
//   		std::cout << "Nearest Samples to given query are:" <<endl;
//   		for (auto it =ranks.rbegin(); it !=ranks.rend(); it++) { // calls a_map.begin() and a_map.end()
// 			if(k>0){
// 				std::cout <<k<<" Similarity "<<it->first <<" Sample column " <<it->second << '\n';	
// 				k--;
// 			}
    		
//   		}
// 	}
// 	end = get_usecs();
// 	dur = ((double)(end-start))/1000000;
// 	printf("\nQuery Search Time = %f\n",dur);
// 	long tend = get_usecs();
// 	double tdur = ((double)(tend-tstart))/1000000;
// 	printf("Total Program Time = %f\n",tdur);
//   	return 0;
// }
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
 //WITH ALGO!!!**

#include <random>
#include <iostream>
#include <string.h>
#include <array>
#include <unordered_map>
#include <string>
#include <stdlib.h>     
#include <time.h>   
#include <vector>
#include "lsh.h"
#include <sys/time.h>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <map>

using namespace std;
// to do add it to the lsh header file Hamming distance help taken from stackoverflow.
int GetHammingDistance(const std::string a, const std::string b)
{
    std::string::const_iterator a_it = a.begin();
    std::string::const_iterator b_it = b.begin();
    std::string::const_iterator a_end = a.end();
    std::string::const_iterator b_end = b.end();
    int value = 0;
    while (a_it != a_end && b_it != b_end)
    {
        if (*a_it != *b_it){ 
        	++value;
        }
        ++a_it;
        ++b_it;
    }
    return value;
}
//=====================================================================================
//=====================================================================================

// here we are using the concept that for it to be similar,! mean - U should be small

// RETURNS 1 for good hyp plane and 0 for bad.
int mann_Whitney_U_test(vector<int> near,vector<int> far)
{
	int near_size = far.size();
	int far_size = near.size();

	int tot_size = near_size + far_size;
	vector<pair <int,int> > arr;
	// two arrays - combined into one and sorted into ascending order
	// keep information about where did it come from
	for (int i = 0; i < near_size; ++i)
		arr.push_back(make_pair(near[i],2));
	for (int i = 0; i < far_size; ++i)
		arr.push_back(make_pair(far[i],1));
	sort(arr.begin(),arr.end());
	// ranks of first sample elements are summarized

	int i = 0,count = 0,num = arr[0].first;
	vector<int> rank;

	// after sorting, each element replaced by its rank 1 - N+M
	while(i < tot_size)
	{
		if(arr[i].first == num)
		{
			count++;
		}
		if(arr[i].first != num)
		{
			int avg = ((2*i)-count+1)/2.0;
			for (int k = 0; k < count; ++k)
				rank.push_back(avg);

			count = 1;
			num = arr[i].first;
		}
		i++;
	}
	int sum1 = 0, sum2 = 0;
	for (int i = 0; i < tot_size; ++i)
	{
		if(arr[i].second == 1)
			sum1 += rank[i];

		else if(arr[i].second == 2)
			sum2 += rank[i];
	}
	// pattern - critical values - 8 8 7 8 8 7... so i subtract 7.66

	int U1, U2, U;
	U1 = (near_size*far_size) + ((near_size*(near_size+1))/2) - sum1;
	U2 = (near_size*far_size)+ ((far_size*(far_size+1))/2) - sum2;

	U = min(U1,U2);
	int mean_U = 0.5 * far_size * near_size;
	//float alpha = 0.05; std deviation
	// U calculated U = N*M + (N(N+1))/2 - summ(rank(xi))
	float std_dev = ((far_size * near_size)*(near_size + far_size + 1))/12.0;
	
	// U1+U2 = M*N
	// if( u close to mean_U)
	//		medians of X&Y are close
	if((mean_U - U) <= (0.5*std_dev))
	{
		//its a bad hyperplane
		return 0;
	}
	else{
		//its gooooddd
		return 1;
	}
}

//=====================================================================================
//=====================================================================================
//=====================================================================================

// MANN WHITNEY U TEST 2 !! Ismei make fixed size arrays. and then use the critical value to check.
int mann_Whitney_U_test2(vector<int> near,vector<int> far)
{
	int near_size = far.size();
	int far_size = near.size();

	int tot_size = near_size + far_size;
	vector<pair <int,int> > arr;
	// two arrays - combined into one and sorted into ascending order
	// keep information about where did it come from
	for (int i = 0; i < near_size; ++i)
		arr.push_back(make_pair(near[i],2));
	for (int i = 0; i < far_size; ++i)
		arr.push_back(make_pair(far[i],1));
	sort(arr.begin(),arr.end());
	// ranks of first sample elements are summarized

	int i = 0,count = 0,num = arr[0].first;
	vector<int> rank;

	// after sorting, each element replaced by its rank 1 - N+M
	while(i < tot_size)
	{
		if(arr[i].first == num)
		{
			count++;
		}
		if(arr[i].first != num)
		{
			int avg = ((2*i)-count+1)/2.0;
			for (int k = 0; k < count; ++k)
				rank.push_back(avg);

			count = 1;
			num = arr[i].first;
		}
		i++;
	}
	int sum1 = 0, sum2 = 0;
	for (int i = 0; i < tot_size; ++i)
	{
		if(arr[i].second == 1)
			sum1 += rank[i];

		else if(arr[i].second == 2)
			sum2 += rank[i];
	}

	// U calculated U = N*M + (N(N+1))/2 - summ(rank(xi))
	int U1, U2, U;
	U1 = (near_size*far_size) + ((near_size*(near_size+1))/2) - sum1;
	U2 = (near_size*far_size)+ ((far_size*(far_size+1))/2) - sum2;

	U = min(U1,U2);
	//int mean_U = 0.5*far_size*near_size;
	// if( u close to mean_U)
	//		medians of X&Y are close
	//int alpha = 0.05;
	
	// pattern - critical values - 8 8 7 8 8 7... so i subtract 7.66
	
	return 0;
}

//=====================================================================================
//=====================================================================================
int good_bad_Hyperplane(float **dataset,float **matrix,int p,int rows,int col,int num_hplane)
{
	int num_pts  = 500;
	//int point = 0,count = 0;
	int point1 = 0;
	int point2 = 0;
	srand (time(NULL));
	float sum = 0, hamm_dist=0;
	
	vector< pair <float,int> > array;
	
	for (int i = 0; i < num_pts; ++i)
	{
		point1 = rand()%col;
		point2 = rand()%col;
		sum = 0,hamm_dist=0;
		for (int j = 0; j < col; ++j)
		{
			// sum == eucledian distance between pt1 and pt0.
			// need to try with cosine distance
			sum += (dataset[j][point1] - dataset[j][point2])*(dataset[j][point1] - dataset[j][point2]);
		}
		
		// hamm dist - should be zero for closer points and 1 for far points.. either 0 or 1
		hamm_dist = abs(matrix[num_hplane][point1] - matrix[num_hplane][point2]);
		array.push_back( make_pair(sum,hamm_dist) );
	}
	sort(array.begin(),array.end());
	//SORTED
	vector<int> near;
	vector<int> far;
	for (int i = 0; i < num_pts; ++i)
	{
		if(array[i].second == 0)
		{
			near.push_back(array[i].first);
		}
		else
		{
			far.push_back(array[i].first);
		}
	}
	// WILKOXYN TEST: If cluster is cut - We'll get high P-values
	// perform rank whitney U test
	
	//their hashes compared..
	return mann_Whitney_U_test(near,far);	
}

long get_usecs(void)
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
	long tstart = get_usecs();
	
	//VARIABLE STUFF
	string path = "../../../../../akarsha.txt";
	int p = 50, rows = 24706,col = 50000;
	string path1 = "../../data/query_738.txt";
	int k = 50;
	

	printf("p = %d\n",p);
	printf("Dataset used:  %s\n",path.c_str());
	long start = get_usecs();
	float ** dataset = lsh::read_data(path);
	long end = get_usecs();
	double dur = ((double)(end-start))/1000000;
	printf("Dataset Load Time = %f\n",dur);
	//  TODO no rows are hardcorded
	// p = number of hyperplanes created ie. k !!
	// col = dimensionality
	// rows = number of datapoints 
	
	float ** trans_data;
	trans_data = new float*[col];
	for(int i = 0; i < col; ++i)
    	trans_data[i] = new float[rows];
	//==========HYPERPLANE========================== 
	for(int i = 0; i < rows; ++i)
		for(int j = 0; j < col; ++j)
            trans_data[j][i]=dataset[i][j];
	
	start = get_usecs();
	float ** hyperplane= lsh::hyper_plane(p,rows);    
	end = get_usecs();
	dur = ((double)(end-start))/1000000;
	printf(" Hyperplane  Time = %f\n",dur);
	
	//=======HASH MATRIX================================
	start = get_usecs();
	float ** hash_matrix = lsh::hash_matrix(dataset,hyperplane,p,col,rows);
	end = get_usecs();
	dur = ((double)(end-start))/1000000;
    printf("Hash matrix  Time = %f\n",dur);
   
    //MY TESTING!!
    int test_matrix[p];
	int hyp_final=0;
	for (int i = 0; i < p; ++i)
	{
		test_matrix[i] = good_bad_Hyperplane(dataset,hash_matrix,p,rows,col,i);
		if(test_matrix[i] == 1)
			hyp_final ++;
	}

	//========HASH TABLES================================
	start = get_usecs();
	std::unordered_map<std::string, std::vector<int> > hash_tables;
	hash_tables = lsh::hash_table(hash_matrix,test_matrix,rows,col,p); 
	end = get_usecs();
	dur = ((double)(end-start))/1000000;
	printf("Hash table  Time = %f\n",dur);
	
	//=========QUERY===============================
	std::unordered_map<std::string, std::vector<int> >:: iterator itr;
	start = get_usecs();

	printf("QueryDataset: %s\n",path1.c_str());
    float query_trans[1][rows];
    float ** query=lsh::read_data(path1);	
    for (int i = 0; i < rows; i++)
    {
		query_trans[0][i]=query[i][0];
    }

	// generating the hash code for the query
 	std::string code="";
	hash_matrix = lsh::hash_matrix(query,hyperplane,p,1,rows);
	for (int i = 0; i < p; i++)
	{
		if(test_matrix[i] == 1){ //threshold = 10
			code.append(std::to_string((int)hash_matrix[i][0]));
		}
	}
	// Doing a K nearest search based on string compare function TODO cosine etc
	int rank[hash_tables.size()];
	cout << "hash_table size: "<< hash_tables.size() << " | "<<"hyp_final = "<<hyp_final << endl;
	int samples[hash_tables.size()];
	std::map<float, int> ranks;
	std::string code1 = code;
	// std::cout<<cosine_similarity(trans_data[0],query_trans[0],rows)<<endl;
	int j = 0;
	if (hash_tables.find(code) == hash_tables.end())
	{
		// cout << "No exact match";
		int i = 0;
		for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
		{
			rank[i] = GetHammingDistance(code1,itr->first);
			i++;
		}
		i=0;
		for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
		{
			if(rank[i] < 10)
			{
				// cout <<itr->first;
				std::vector<int> val = itr->second;
			    while (!val.empty()) 
			    {
		   		    auto it = val.rbegin();
		    	    // printf("%d\n", *it);
		    		samples[j]=*it;
		    		val.pop_back();
		    		j++;
		  		}
			}
			i++;
		}
		printf("ehaa\n");
		for (int k = 0; k < j; k++)
		{
			ranks[cosine_similarity(trans_data[samples[k]],query_trans[0],rows)]=samples[k];
			cout << samples[k]<<"\t"<<cosine_similarity(trans_data[k],query_trans[0],rows)<<"\n";
		}
  		std::cout << "Nearest Samples to given query are:" <<endl;
  		for (auto it =ranks.rbegin(); it !=ranks.rend(); it++) { // calls a_map.begin() and a_map.end()
			if(k > 0 ){
				std::cout <<"Similarity "<<it->first <<" Sample column " <<it->second << '\n';	
				k--;
			}
    		
  		}
	}
	// Key found part
	else
	{
		printf("lalaa\n");
		int temp = k;
		std::vector<int> val = hash_tables[code1];
	    while (!val.empty()) {
   			 auto it = val.rbegin();
    		 // printf("%d\n", *it);
    		 samples[j]=*it;
    		 val.pop_back();
    		 j++;
  		}
  		if(j == temp)
  		{
  			for (int k = 0; k < j; k++)
			{
				ranks[cosine_similarity(trans_data[samples[k]],query_trans[0],rows)]=samples[k];
			}
  		}
  		else
  		{
			int i=0;
			for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
			{
				if(itr->first == code1)
					continue;
				else
				{
					rank[i] = GetHammingDistance(code1,itr->first);
					i++;
				}
			}
			i = 0;
			for (itr = hash_tables.begin(); itr != hash_tables.end(); itr++)
			{
				if(rank[i] < 10)
				{
					// cout <<itr->first;
					std::vector<int> val = itr->second;
				    while (!val.empty()) {
			   			 auto it = val.rbegin();
			    		 // printf("%d\n", *it);
			    		 samples[j] =* it;
			    		 val.pop_back();
			    		 j++;
			  		}
				}
				i++;
			}
			for (int k = 0; k < j; k++)
			{
				ranks[cosine_similarity(trans_data[samples[k]],query_trans[0],rows)]=samples[k];
				//cout << samples[k]<<"\t"<<cosine_similarity(trans_data[k],query_trans[0],rows)<<"\n";
			}

  		}
  		std::cout << "Nearest Samples to given query are:" <<endl;
  		for (auto it =ranks.rbegin(); it !=ranks.rend(); it++) { // calls a_map.begin() and a_map.end()
			if(k>0){
				cout << "k = "<< k ;
				std::cout <<"  Similarity "<< it->first <<" Sample column " <<it->second << '\n';	
				k--;
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
