/**
 *	lsh_cpu.H
 *
 */
// #ifndef _lsh_H_INCLUDED
// #define _lsh_H_INCLUDED
// #include <string.h>
// #include <iostream>
// #include <unordered_map>
// #include <array>
// #include <vector>
// #include <string>
// namespace lsh{
// 	/**
// 	 *	@param int text - Integer Text to print on the screen.
// 	 */
// 	void example_printf(int text);
// 	/**
// 	 *	@param text - File path to read the file and return the matrix data
// 	 */
// 	float ** read_data(std::string path);
// 	/**
// 	 *	@param text - File path to read the file and return the matrix data
// 	 */
// 	float ** hyper_plane(int p, int rows);	
// 	/**
// 	 *	@param text - Multiply hyperplane and data matrix to generate the hash code
// 	 */
// 	float ** hash_matrix(float **data , float ** hyperplane,int p ,int n, int rows);

// 	std::unordered_map<std::string, std::vector<int>> hash_table(float **matrix, int col,int p);

// }

// #endif


/**
 *	lsh_cpu.H
 *	MINE :D 
 */

#ifndef _lsh_H_INCLUDED
#define _lsh_H_INCLUDED
#include <string.h>
#include <iostream>
#include <unordered_map>
#include <array>
#include <vector>
#include <string>
namespace lsh{
	/**
	 *	@param int text - Integer Text to print on the screen.
	 */
	void example_printf(int text);
	/**
	 *	@param text - File path to read the file and return the matrix data
	 */
	float ** read_data(std::string path);

	//	@param text - File path to read the file and return the matrix data
	 
	float ** hyper_plane(int p, int rows);	
	/**
	 *	@param text - Multiply hyperplane and data matrix to generate the hash code
	 */
	// float ** hash_matrix(float **data , float ** hyperplane,int p ,int n, int rows);
	float ** hash_matrix(int whichone,float **dataset, float **hyperplane,int p ,int col, int rows);

	std::unordered_map< std::string, std::vector<int> > hash_table(float **matrix, int* test_matrix, int rows, int col,int p);

}

#endif
