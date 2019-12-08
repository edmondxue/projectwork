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

   Th Th Th Th ... Th
   O  O  O  O  ... P
   O  O  O  O  ... P
   Tx Tx Tx Tx ... Tx

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
	const int nx = (int) (len / h);
	const int ny = (int) (width / h);
		
	A.Resize(ny, nx);

	//set const coeff
	const double coeff = 1 / pow(h, 2.0);

	//form initial product vector b
	this->b.resize((ny - 1) * nx, 0);

	//form initial soln guess vector x
	this->x.resize((ny - 1) * nx, 1);

	//loop through all points to be solved
	for (int i = 0; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			//translation of i,j position to A's i,j

			A.AddEntry(to1D(i, j, nx), to1D(i - 1, j, nx), 1);
			A.AddEntry(to1D(i, j, nx), to1D(i + 1, j, nx), 1);

			//check if lower isothermal boundary
			if (j - 1 > 0)
			{
				//if not, normally add to A
				A.AddEntry(to1D(i, j, nx), to1D(i, j-1, nx), 1);
			}
			else
			{
				//if it is, add Tx to b
				double Tx = -Tc * (exp(-10 * pow(len/nx * i - (len / 2), 2.0)) - 2);
				this->b[to1D(i, j, nx)] -= Tx;
			}

			//check if upper isothermal boundary
			if (j + 1 < ny)
			{
				//if not, normally add to A
				A.AddEntry(to1D(i, j, nx), to1D(i, j + 1, nx), 1);
			}
			else
			{
				//if it is, add Th to b
				this->b[to1D(i, j, nx)] -= Th;

			}
			//always account for term of current point
			A.AddEntry(to1D(i, j, nx), to1D(i, j, nx),  -4);
		}
	}

	return 0;
}

/* Method to solve system using CGsolver */
int HeatEquation2D::Solve(std::string soln_prefix)
{
	//set tolerance and call to CGSolver
	double const tol = 1.0 * pow((double)10.0, -5);

	//implement negative definite (-A)u = -b

	this->A.MulConst(-1.0);
	this->b = constvec_mult(-1.0, this->b);
	
	int iter = CGSolver(this->A, this->b, this->x, tol, soln_prefix, *this);


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
std::vector<double> HeatEquation2D::getTemps() const
{
	std::vector<double> temps = {this->Tc, this->Th};
	
	return temps;
}

/* Method to get length, width, h for the system;*/
std::vector<double> HeatEquation2D::getDims() const
{
	std::vector<double> dims = { this->len, this->width, this->h};

	return dims;
}