#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>


#include "COO2CSR.hpp"
#include "CGSolver.hpp"


/* Loads a matrix from a file specified at the command line, converts the matrix
 * to CSR format, runs the CG solver function with a starting guess of ones for the solution
 * and zeros for the right hand side, and writes the solution to the specified file name
 */

int main(int argc, char *argv[])
{
	// Ensure at least 3 parameters passed
	if (argc < 3)
	{
		std::cout << "Usage: " << std::endl;
		std::cout << "$ ./main matrix2.txt solution2.txt" << std::endl;
		return 0;
	}

	//initialize
	std::string matrixfile, solutionfile;

	std::vector<double> val;
	std::vector<int> i_idx;
	std::vector<int> j_idx;
	int rowsize, colsize;


	matrixfile = argv[1];
	solutionfile = argv[2];

	//read in the matrix file
	std::ifstream f(matrixfile.c_str());
	if(f.is_open())
	{

		int row_ind, col_ind;
		float value;

		//read first row into sizes
		f >> rowsize >> colsize;

		while(f >> row_ind >> col_ind >> value)
		{
			val.push_back(value);
			i_idx.push_back(row_ind);
			j_idx.push_back(col_ind);

		}
		f.close();

	}
	else
	{
		std::cerr << "ERROR: Failed to open file" << std::endl;
		return 0;
	}


	//after reading matrix, call COO2CSR
	//convert matrix format
	//i_idx becomes rowptr, j_idx becomes col
	COO2CSR(val, i_idx, j_idx);

	//solve: tol of 1e-5, x soln guess 1's, b rhs guess 0's
	double const tol = 1.0 * pow((double) 10.0, -5);
	std::vector<double> bvec (rowsize, 0.0);
	std::vector<double> xvec (rowsize, 1.0);


	int iter = CGSolver(val, i_idx, j_idx, bvec, xvec, tol);

	//solution vector now in xvec
	//write it to output file

	std::ofstream outf(solutionfile.c_str());
	if(outf.is_open())
	{
		for (unsigned int i = 0; i < xvec.size(); i++)
		{
			//fix output format
			outf << std::setprecision(5) << xvec[i] << "\n";
		}
		outf.close();
	}
	else
	{
		std::cerr << "ERROR: Failed to open output file" << std::endl;
		return 0;
	}

	//output success or failure message
	if (iter != -1)
	{
		std::cout << "SUCCESS: CG solver converged in " << iter << " iterations."<< std::endl;
	}
	else
	{
		std::cout << "Solver Didn't Converge" << std::endl;
	}

	return 0;
}

