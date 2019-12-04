#ifndef HEAT_HPP
#define HEAT_HPP

#include <string>
#include <vector>

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
	int Setup(std::string inputfile);

    /* Method to solve system using CGsolver */
    int Solve(std::string soln_prefix);

	/* Method to get Tc and Th */
	std::vector<double> getTemps();

    /* TODO: Add any additional public methods you need */

};

#endif /* HEAT_HPP */
