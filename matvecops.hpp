#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>

/*  additional functions to perform common vector and matrix
 *  operations that occur in the CG algorithm.
 */

//returns difference between vectors
std::vector<double> vec_subtract(std::vector<double> const& x, std::vector<double> const& y);

//returns sum between vectors
std::vector<double> vec_add(std::vector<double> const& x, std::vector<double> const& y);

//returns multiplication product between vectors
//transposing one
std::vector<double> matvec_mult(std::vector<double> const& val,
								std::vector<int> const& row_ptr,
								std::vector<int> const& col_idx,
								std::vector<double> const& x);

//return vector product of CSR matrix and constant
void matconst_mult(std::vector<double>& val, std::vector<int> const& row_ptr,
	std::vector<int> const& col_idx, double const conval);

//return dot product of two vectors, transposing one
double dot_prod(std::vector<double> const& x, std::vector<double> const& y);

//return product of constant and a vector
std::vector<double> constvec_mult(double const con, std::vector<double> const& x);

//calculate the l2 norm of a vector
double L2norm(std::vector<double> const& x);


#endif /* MATVECOPS_HPP */
