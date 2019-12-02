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
		//ncols = num unknowns = num interior + num periodic/2
		//nrows = num equations = num unknowns?
		const int nx = len / h;
		const int ny = width / h;
		
		const int n_unk = (nx - 1) * (ny - 1) + (nx - 1);
		A.Resize(n_unk, n_unk);


		////first, resize the Sparse Matrix
		////nrows and ncols are actually ndivs, add 1 to capture indices
		//const int ncols = len / h;
		//const int nrows = width / h;

		//A.Resize(nrows + 1, ncols + 1);


		//two types of BC's: periodic and isothermal

		//first, isothermal:
		
		//top layer is all at Th
		//for (int col_ind = 0; col_ind < len; col_ind++)
		//{
		//	A.AddEntry(0, col_ind, Th);
		//}

		////bottom layer defined as function Tx
		//double Tx;

		//for (int col_ind = 0; col_ind < len; col_ind++)
		//{
		//	
		//}
		
		//second, periodic:


		//set const coeff
		const double coeff = 1 / pow(h, 2.0);
		////loop through all interior points i,j
		//for (int i = 0; i < nrows + 1; i++)
		//{
		//	for (int j = 0; j < ncols + 1; j++)
		//	{
		//		//for isothermal BC's
		//		//top
		//		if (i == 0)
		//		{
		//			//at the j'th 
		//			A.AddEntry(j, j, Th);
		//		}
		//		//bottom
		//		if (i == nrows)
		//		{
		//			double Tx = -Tc * (exp(-10 * pow(j - (len / 2), 2.0)) - 2);
		//			A.AddEntry(width - 1, j, Tx);
		//		}


		//		//check boundaries
		//		if((i-1) < 0
		//		A.AddEntry(coeff, )
		//	}
		//}

		//form initial product vector b
		std::fill(b.begin(), b.end(), 0);

		//form initial soln guess vector x
		std::fill(x.begin(), x.end(), 1);


		for (int j = 1; j < ny; j++)
		{
			for (int i = 1; i < nx; i++)
			{
				// if i - 1 = 0, i - 1 periodic BC
				// if i + 1 = nx, i + 1 periodic BC
				if (i)
				{
					A.AddEntry();
				}
				//if j - 1 = 0, j - 1 bottom isotherm BC
				if ((j - 1) == 0)
				{
					double Tx = -Tc * (exp(-10 * pow(j - (len / 2), 2.0)) - 2);
					b[nx * j + i] -= coeff * Tx;
				}
				//if j + 1 = ny, j + 1 top isotherm BC
				if ((j + 1) == ny)
				{
					b[nx * j + i] -= coeff * Th;
				}
				

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