#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include <vector>

#include "heat.hpp"
#include "sparse.hpp"

/* Function that implements the CG algorithm for a linear system
 *
 * Ax = b
 *
 * where A is in CSR format.  The starting guess for the solution
 * is provided in x, and the solver runs a maximum number of iterations
 * equal to the size of the linear system.  Function returns the
 * number of iterations to converge the solution to the specified
 * tolerance, or -1 if the solver did not converge.
 */
int CGSolver(SparseMatrix& mat, std::vector<double> const& b, std::vector<double>& x,
	const double tol, const std::string soln_prefix, HeatEquation2D const & sys);

/*Function that fills in boundary conditions left out of A and prints the 
solution to a new file*/
int printSolnFile(const std::string soln_prefix, std::vector<double> const& x,
	const int niter, SparseMatrix& mat, HeatEquation2D const & sys);

#endif /* CGSOLVER_HPP */
