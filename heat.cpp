#include <cmath>
#include <fstream>
#include <iostream>
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
	double len, width, h;
	double Tc, Th;


	/* TODO: Add any additional private data attributes and/or methods you need */

public:
	/* Method to setup Ax=b system */
	int Setup(std::string inputfile)
	{
		std::ifstream f(inputfile.c_str());
		if (f.is_open())
		{

			//read first row into geometry
			f >> this->len >> this->width >> this->h;

			//read second row into temps
			f >> this->Tc >> this->Th;

			//close file
			f.close();

		}
		else
		{
			std::cerr << "ERROR: Failed to open file" << std::endl;
			return 0;
		}
		
		//first, resize the Sparse Matrix
		//nrows and ncols are actually ndivs, add 1 to capture indices
		const int nrows = len / h;
		const int ncols = width / h;
		
		A.Resize(nrows, ncols);


		//two types of BC's: periodic and isothermal

		//first, isothermal:
		
		//top layer is all at Th
		for (int col_ind = 0; col_ind < len; col_ind++)
		{
			A.AddEntry(0, col_ind, Th);
		}

		//bottom layer defined as function Tx
		double Tx;

		for (int col_ind = 0; col_ind < len; col_ind++)
		{
			Tx = -Tc * (exp(-10 * pow(col_ind - (len / 2), 2.0)) - 2);
			A.AddEntry(width-1, col_ind, Tx);
		}
		
		//second, periodic:


		//set const coeff
		const 
		//loop through all points i,j
		for (int i = 0; i < ncols + 1; i++)
		{
			for (int j = 0; j < nrows + 1; j++)
			{
				A.AddEntry()
			}
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
	x
	/* TODO: Add any additional public methods you need */

};