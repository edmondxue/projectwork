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
	//nrows and ncols are actually ndivs, add 1 to capture indices
	//ncols = num unknowns = num interior + num periodic/2
	//nrows = num equations = num unknowns?
	const int nx = (int) (len / h);
	const int ny = (int) (width / h);
		
	//const int n_unk = (int) ((nx - 1) * (ny - 1) + (nx - 1));
	A.Resize(ny, nx);



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
	this->b.resize((ny - 1) * nx, 0);

	std::cout << (ny - 1) * nx << " " << b.size();
	//std::fill(this->b.begin(), this->b.end(), 0);


	//form initial soln guess vector x
	this->x.resize((ny - 1) * nx, 1);
	//std::fill(this->x.begin(), this->x.end(), 1);

	//vars for adding to A
	int Ai, Aj; // REMOVE

	std::cout << A.getDims()[0] << " " << A.getDims()[1] << "\n";
	std::cout << nx << " " << ny << "\n";

	////loop through all points to be solved
	//for (int i = 0; i < nx; i++)
	//{
	//	for (int j = 1; j < ny; j++)
	//	{
	//		//translation of i,j position to A's i,j
	//		Ai = nx * (j - 1) + i;
	//		Aj = Ai;


	//		// if i = 0, curr on top of periodic BC
	//		if (i == 0)
	//		{
	//			//right side normal
	//			A.AddEntry(Ai, Aj + 1, 1);
	//			//left side loops over to i = nx -1
	//			A.AddEntry(Ai, Aj + (nx - 1), 1);
	//		}
	//		// if i = nx - 1, right refers to periodic BC
	//		else if (i == nx - 1)
	//		{
	//			//left side normal
	//			A.AddEntry(Ai, Aj - 1, 1);
	//			//right side loops back to i = 0
	//			A.AddEntry(Ai, Aj - (nx-1), 1);
	//		}
	//		//otherwise normal left/right
	//		else
	//		{
	//			A.AddEntry(Ai, Aj - 1, 1);
	//			A.AddEntry(Ai, Aj + 1, 1);
	//		}


	//		//normal top/bottom,  not near isotherm
	//		if  ((j - 1) != 0 and (j + 1) != ny)
	//		{
	//			A.AddEntry(Ai, Aj + nx, 1);
	//			A.AddEntry(Ai, Aj - nx, 1);
	//		}
	//		//if j - 1 = 0, j - 1 bottom pt isotherm BC
	//		else if ((j - 1) == 0)
	//		{
	//			//add Tx to b, don't add any to A
	//			double Tx = -Tc * (exp(-10 * pow(i - (len / 2), 2.0)) - 2);
	//			//b[Ai] -= coeff * Tx;
	//			b[Ai] -= Tx;
	//			//normal top
	//			A.AddEntry(Ai, Aj + nx, 1);
	//		}
	//		//if j + 1 = ny, j + 1 top pt isotherm BC
	//		else if ((j + 1) == ny)
	//		{
	//			//add Th to b, don't add any to A
	//			//b[Ai] -= coeff * Th;
	//			b[Ai] -=  Th;
	//			//normal bottom
	//			A.AddEntry(Ai, Aj - nx, 1);
	//		}
	//				
	//		//std::cout << "vdf3";
	//		//always account for term of current point
	//		A.AddEntry(Ai, Aj, -4*1);


	//	}
	//}


	//loop through all points to be solved
	for (int i = 0; i < nx; i++)
	{
		for (int j = 1; j < ny; j++)
		{
			//translation of i,j position to A's i,j
			Ai = nx * (j - 1) + i;
			Aj = Ai;

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
				double Tx = -Tc * (exp(-10 * pow(i - (len / 2), 2.0)) - 2);
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
				this->b[Ai] -= Th;

			}
			//always account for term of current point
			A.AddEntry(to1D(i, j, nx), to1D(i, j, nx),  -4);


		}
	}

	//A.printIt();
	//
	//std::cout <<  "\n ";

	//for (unsigned int i = 0; i < b.size(); i++)
	//{
	//	std::cout << b[i] << " ";
	//}

	//for (unsigned int i = 0; i < x.size(); i++)
	//{
	//	std::cout << x[i] << " ";
	//}

	


	return 0;
}

/* Method to solve system using CGsolver */
int HeatEquation2D::Solve(std::string soln_prefix)
{

	//set tolerance and call to CGSolver
	double const tol = 1.0 * pow((double)10.0, -5);
	//int CGSolver(SparseMatrix& mat, std::vector<double> const& b, std::vector<double>& x, 
	//const double tol, const std::string soln_prefix, HeatEquation2D& sys)

	//implement negative definite (-A)u = -b

	this->A.MulConst(-1.0);
	this->b = constvec_mult(-1.0, this->b);


	/*std::cout << "\n\n";
	A.printIt();

	std::cout << "\n";

	for (unsigned int i = 0; i < b.size(); i++)
	{
		std::cout << b[i] << " ";
	}
	std::cout << "\n";

	for (unsigned int i = 0; i < x.size(); i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << "\n";*/
	
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




	/* TODO: Add any additional public methods you need */
