
#include <iostream>
#include <string.h>
#include <array>
#include <unordered_map>
#include <string>
#include <vector>
#include <sys/time.h>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <map>
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

using namespace std;

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
	string path = "data/dataset_2048.txt";
	string line;
	int row=0,col=0;
	float** data;
	float** qdata;	
	data = new float*[60000];
	for(int i = 0; i < 60000; ++i)
    	data[i] = new float[2049];
	
	ifstream myfile (path);
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
	row=58307,col=2049;
	float ** trans_data;
	trans_data = new float*[2049];
	for(int i = 0; i < 2049; ++i)
    	trans_data[i] = new float[58307];
	
	for(int i = 0; i < row; ++i)
	{
		for(int j = 0; j < col; ++j)
        {
            trans_data[j][i]=data[i][j];
     	}
	}
	string path1 = "data/query_738.txt";
    float query_trans[1][row];
    qdata = new float*[60000];
	for(int i = 0; i < 60000; ++i)
    	qdata[i] = new float[2049];
	
	row=0,col=0;
	ifstream myfileq (path1);
	if (myfileq.is_open())
	{
		while ( getline (myfileq,line) )
		{
			std::istringstream iss(line);
			int val;
			col=0;
			while ( iss >> val)
			{
			
				qdata[row][col]=val;
				col++;
			}
	  		row++;
		}
		myfileq.close();
	}
	else cout << "Unable to open file"; 

	row=58307,col=2049;
    for (int i = 0; i < row; i++)
    {
    	for (int j = 0; j < 1; j++)
    	{
    		query_trans[j][i]=qdata[i][j];
    	}
     	 
    }
    int x = 40;
    std::map<float, int> ranks;
    for (int k = 0; k < col; k++)
	{
		ranks[cosine_similarity(trans_data[k],query_trans[0],row)]=k;
	}
	for (auto it =ranks.rbegin(); it !=ranks.rend(); it++) 
	{ 		// calls a_map.begin() and a_map.end()
		std::cout << x<<" Similarity "<<it->first <<" Sample column " <<it->second << '\n';
		x--;
  	}
    
	return 0;
}