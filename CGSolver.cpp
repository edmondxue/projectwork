#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "CGSolver.hpp"
#include "heat.hpp"
#include "matvecops.hpp"
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
			const double tol, const std::string soln_prefix, HeatEquation2D const& sys)
{
	//initialize
	std::vector<double> u, u_new, r, r_new, p, p_new;
	double L2normr0, L2normr, alpha, beta;
	bool converged = false;

	//begin CG algo
	u = x;
	r = vec_subtract(b, mat.MulVec(u));
	L2normr0 = L2norm(r);
	p = r;

	int niter = 0;
	const int nitermax = (int) b.size();

	while(niter < nitermax)
	{
		niter += 1;
		alpha = dot_prod(r, r)/dot_prod(p, mat.MulVec(p));
		u_new = vec_add(u, constvec_mult(alpha,p));
		r_new = vec_subtract(r, constvec_mult(alpha, mat.MulVec(p)));
		L2normr = L2norm(r_new);

		if(L2normr/L2normr0 < tol)
		{
			converged = true;
			break;
		}

		beta = dot_prod(r_new, r_new)/dot_prod(r,r);
		p_new = vec_add(r_new, constvec_mult(beta,p));

		//update variables
		r = r_new;
		p = p_new;
		u = u_new;
		x = u_new;

		//on every 10th iteration, print out solution file
		if (niter % 10 == 0)
		{
			printSolnFile(soln_prefix, x, niter, sys);
		}


	}

	//print the last iteration
	printSolnFile(soln_prefix, x, niter);


	//if converges, return number of iterations
	if (converged)
	{
		return niter;
	}
	else
	{
		//if doesn't converge, return -1
		return -1;
	}

	
}

int printSolnFile(const std::string soln_prefix,
	std::vector<double> const& x, const int niter, HeatEquation2D const& sys)
{
	//format the niter for filename
	std::stringstream niter_str;
	niter_str << std::setw(4) << std::setfill('0') << niter;

	//use the x vector to create a solution that contains the isothermal boundary points
	double Tc = sys.getTemps()[0];
	double Th = sys.getTemps()[1];

	for (int x = 0; x < sys.)
	double Tx = -Tc * (exp(-10 * pow(j - (len / 2), 2.0)) - 2);

	std::ofstream outf(soln_prefix.c_str() + niter_str.str());
	if (outf.is_open())
	{
		for (unsigned int i = 0; i < x.size(); i++)
		{
			//fix output format
			outf << std::setprecision(5) << x[i] << "\n";
		}
		outf.close();
	}
	else
	{
		std::cerr << "ERROR: Failed to open output file" << std::endl;
		return 1;
	}
	return 0;
}
