#ifndef HEAT_HPP
#define HEAT_HPP

#include <string>
#include <vector>

#include "sparse.hpp"

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
	std::vector<double> getTemps() const;

	/* Method to get length, width, h for the system;*/
	std::vector<double> getDims() const;



    /* TODO: Add any additional public methods you need */

};

#endif /* HEAT_HPP */
