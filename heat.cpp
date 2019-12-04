#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "heat.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"


/* Method to setup Ax=b system */

/* A is built to remove calculations for those already known points. This includes
   all isothermal boundary condition points and 1/2 of the periodic boundary
   condition points.

   Example:

   Th Th Th Th Th
   P  O  O  O  O
   P  O  O  O  O
   Tx Tx Tx Tx Tx

   where those points marked O are included in A to be calculated, and other markings 
   are left out to be calculated later.
*/

int HeatEquation2D::Setup(std::string inputfile)
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
		return 1;
	}
		
	//first, resize the Sparse Matrix
	//nrows and ncols are actually ndivs, add 1 to capture indices
	//ncols = num unknowns = num interior + num periodic/2
	//nrows = num equations = num unknowns?
	const int nx = (int) (len / h);
	const int ny = (int) (width / h);
		
	const int n_unk = (int) ((nx - 1) * (ny - 1) + (nx - 1));
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

	//vars for adding to A
	//int Ai_R, Aj_R, Ai_L, Aj_L, Ai_U, Aj_U, Ai_D, Aj_D;
	int Ai, Aj;

	for (int j = 1; j < ny; j++)
	{
		for (int i = 1; i < nx + 1; i++)
		{
			//translation of i,j position to A's i,j
			Ai = nx * (j - 1) + i;
			Aj = Ai;

			//// if i - 1 = 0, i - 1 periodic BC
			//if ((i - 1) == 0)
			//{
			//	//left side, loop to right side unkn
			//	A.AddEntry();
			//}

			// if i = nx, curr on top of periodic BC
			if (i == nx)
			{
				//right side ref loops to i = 1
				A.AddEntry(Ai, Aj - (nx - 1), coeff);
				//left normal
				A.AddEntry(Ai, Aj - 1, coeff);
			}
			//otherwise normal left/right
			else
			{
				A.AddEntry(Ai, Aj - 1, coeff);
				A.AddEntry(Ai, Aj + 1, coeff);
			}


			//normal top/bottom,  not near isotherm
			if  ((j - 1) != 0 and (j + 1) != ny)
			{
				A.AddEntry(Ai, Aj + nx, coeff);
				A.AddEntry(Ai, Aj - nx, coeff);
			}
			//if j - 1 = 0, j - 1 bottom pt isotherm BC
			else if ((j - 1) == 0)
			{
				//add Tx to b, don't add any to A
				double Tx = -Tc * (exp(-10 * pow(j - (len / 2), 2.0)) - 2);
				b[Ai] -= coeff * Tx;
				//normal top
				A.AddEntry(Ai, Aj + nx, coeff);
			}
			//if j + 1 = ny, j + 1 top pt isotherm BC
			else if ((j + 1) == ny)
			{
				//add Th to b, don't add any to A
				b[Ai] -= coeff * Th;
				//normal bottom
				A.AddEntry(Ai, Aj - nx, coeff);
			}
					
			//always account for term of current point
			A.AddEntry(Ai, Aj, -4*coeff);

		}
	}

	return 0;
}

/* Method to solve system using CGsolver */
int HeatEquation2D::Solve(std::string soln_prefix)
{
	//set tolerance and call to CGSolver
	double const tol = 1.0 * pow((double)10.0, -5);
	//int CGSolver(SparseMatrix mat, std::vector<double> const& b, std::vector<double> & x, const double tol)
	//implement negative definite (-A)u = -b
	int iter = CGSolver(this->A, this->b, this->x, tol, soln_prefix, this);


	//output success or failure message
	if (iter != -1)
	{
		std::cout << "SUCCESS: CG solver converged in " << iter << " iterations." << std::endl;
	}
	else
	{
		std::cout << "Solver Didn't Converge" << std::endl;
	}
	return 0;
}

/* Method to get Tc and Th */
std::vector<double> HeatEquation2D::getTemps()
{
	std::vector<double> temps = {this->Tc, this->Th};
	
	return temps;
}

/* Method to get length, width, h for the system;*/
std::vector<double> HeatEquation2D::getDims()
{
	std::vector<double> dims = { this->len, this->width, this->h};

	return dims;
}

	/* TODO: Add any additional public methods you need */
