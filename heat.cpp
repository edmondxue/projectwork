#include <iostream>
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

		//two types of BC's: periodic and isothermal

		//first isothermal:


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