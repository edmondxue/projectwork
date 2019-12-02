#include <fstream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "heat.hpp"
#include "sparse.hpp"


class HeatEquation2D
{
private:
	SparseMatrix A;
	std::vector<double> b, x;

	/* TODO: Add any additional private data attributes and/or methods you need */

public:
	/* Method to setup Ax=b system */
	int Setup(std::string inputfile)
	{
		std::ifstream f(inputfile.c_str());
		if (f.is_open())
		{

			int row_ind, col_ind;
			float value;

			//read first row into sizes
			f >> rowsize >> colsize;

			while (f >> row_ind >> col_ind >> value)
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
	}

	/* Method to solve system using CGsolver */
	int Solve(std::string soln_prefix)
	{
		//set tolerance and call to CGSolver
		double const tol = 1.0 * pow((double)10.0, -5);
		CGSolver(this->A, tol);
		return 0;
	}

	/* TODO: Add any additional public methods you need */

};