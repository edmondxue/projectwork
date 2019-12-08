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
	const double tol, const std::string soln_prefix, HeatEquation2D const & sys)
{
	//initialize
	std::vector<double> u, u_new, r, r_new, p, p_new;
	double L2normr0, L2normr, alpha, beta;
	bool converged = false;

	//first convert to CSR format
	mat.ConvertToCSR();



	//begin CG algo
	u = x;

	//a
	std::cout << "DEBUG\n\n";

	std::cout << "RESULTS OF  MULVEC \n";
	std::vector <double> a = mat.MulVec(u);

	std::cout << "SIZE OF A: " << a.size() << " \n";
	//for (unsigned int i = 0; i < a.size(); i++)
	//{
	//	std::cout << a[i] << " ";
	//}
	std::cout << "END OF  MULVEC \n";

	//a

	r = vec_subtract(b, mat.MulVec(u));
	L2normr0 = L2norm(r);
	p = r;

	int niter = 0;
	const int nitermax = (int) b.size();

	while(niter < nitermax)
	{
		std::cout << "block 1\n";
		alpha = dot_prod(r, r)/dot_prod(p, mat.MulVec(p));
		u_new = vec_add(u, constvec_mult(alpha,p));
		r_new = vec_subtract(r, constvec_mult(alpha, mat.MulVec(p)));
		L2normr = L2norm(r_new);

		if(L2normr/L2normr0 < tol)
		{
			converged = true;
			break;
		}

		std::cout << "block 2\n";
		std::cout << "set beta\n";
		beta = dot_prod(r_new, r_new)/dot_prod(r,r);

		std::cout << "set p new\n";
		p_new = vec_add(r_new, constvec_mult(beta,p));
		std::cout << "end block 2\n";

		//update variables
		r = r_new;
		p = p_new;
		u = u_new;
		x = u_new;

		//on every 10th iteration, print out solution file
		if (niter % 10 == 0)
		{
			printSolnFile(soln_prefix, x, niter, mat, sys);
		}

		std::cout << " inc niter \n";
		niter += 1;
	}

	//print the last iteration
	printSolnFile(soln_prefix, x, niter, mat, sys);


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

/*Function that fills in boundary conditions left out of A and prints the
solution to a new file*/
int printSolnFile(const std::string soln_prefix, std::vector<double> const& x,
	const int niter, SparseMatrix& mat, HeatEquation2D const & sys)
{
	/*use the x vector to create a solution that 
	contains the isothermal, periodic boundary points*/

	//assign system, matrix properties
	double Tc = sys.getTemps()[0];
	double Th = sys.getTemps()[1];
	double Tx;

	double len = sys.getDims()[0];

	int nrows = mat.getDims()[0];
	int ncols = mat.getDims()[1];

	//std::vector <double> soln((ncols+1) * (nrows+1));
	std::vector <double> soln;

	/*x vector of form: [u(0,1), u(1,1)..]
	what this means is j is 1 greater than its reference

	also, j and i order of looping matters*/

	std::cout << "printing loop\n";
	

	//std::cout << "This is x: \n\n";
	//for (unsigned int i = 0; i < x.size(); i++)
	//{
	//	std::cout << x[i] << " ";
	//}
	//std::cout << " \n\n";
	

	std::cout << nrows << " " << ncols;
	for (int j = 0; j <= nrows; j++)
	{			
		for (int i = 0; i <= ncols; i++)
		{

			//if j = 0, fill in bottom BC Tx's
			if (j == 0)
			{
				Tx = -Tc * (exp(-10 * pow(i - (len / 2), 2.0)) - 2);
				soln.push_back(Tx);
				//soln[to1D(i, j, ncols)] = Tx;
			}
			//if j = nrows, fill in top BC Th's 
			else if (j == nrows)
			{
				soln.push_back(Th);
				//soln[to1D(i, j, ncols)] = Th;
			}
			//otherwise it's an interior row
			else
			{
				////last col interior, copy periodic BC
				//if (i == ncols)
				//{
				//	soln.push_back(x[(j-1) * (ncols) + 1]);
				//}
				//else
				//{
				//	soln.push_back(x[ncols * (j - 1) + i]);
				//}
				soln.push_back(x[to1D(i,j,ncols)]);
				//soln[to1D(i, j, ncols)] = x[to1D(i, j, ncols)];
			}
		}
		//std::cout << " \n\n";
		//for (unsigned int i = 0; i < soln.size(); i++)
		//{
		//	std::cout << soln[i] << " ";
		//}
	}

	std::cout << "finished writing soln\n";
	
	//format the niter for filename
	std::stringstream niter_str;
	niter_str << std::setw(4) << std::setfill('0') << niter;

	std::ofstream outf(soln_prefix.c_str() + niter_str.str() + "neg.txt");
	if (outf.is_open())
	{
		for (unsigned int i = 0; i < soln.size(); i++)
		{
			//fix output format
			outf << std::setprecision(5) << soln[i] << "\n";
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
