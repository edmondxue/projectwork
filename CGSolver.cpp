#include <iostream>
#include <vector>

#include "matvecops.hpp"

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
int CGSolver(std::vector<double> &val,
             std::vector<int>    &row_ptr,
             std::vector<int>    &col_idx,
             std::vector<double> &b,
             std::vector<double> &x,
             double              tol)
{
	//initialize
	std::vector<double> u, u_new, r, r_new, p, p_new;
	double L2normr0, L2normr, alpha, beta;
	bool flag = false;

	//begin CG algo
	u = x;
	r = vec_subtract(b, matvec_mult(val, row_ptr, col_idx, u));
	L2normr0 = L2norm(r);
	p = r;
	int niter = 0;
	const int nitermax = (int) b.size();

	while(niter < nitermax)
	{
		niter += 1;
		alpha = dot_prod(r, r)/dot_prod(p,matvec_mult(val, row_ptr, col_idx, p));
		u_new = vec_add(u, constvec_mult(alpha,p));
		r_new = vec_subtract(r, constvec_mult(alpha, matvec_mult(val, row_ptr, col_idx, p)));
		L2normr = L2norm(r_new);

		if(L2normr/L2normr0 < tol)
		{
			flag = true;
			break;
		}

		beta = dot_prod(r_new, r_new)/dot_prod(r,r);
		p_new = vec_add(r_new, constvec_mult(beta,p));

		//update variables
		r = r_new;
		p = p_new;
		u = u_new;
	}


	//if converges, set solution vector to found solution
	//return number of iterations
	if (flag)
	{
		x = u_new;
		return niter;
	}
	else
	{
		//if doesn't converge, return -1
		return -1;
	}

}