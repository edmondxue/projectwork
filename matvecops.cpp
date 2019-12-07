#include <cmath>
#include <vector>

#include <algorithm>
#include "matvecops.hpp"


//do matrix subtraction between two vectors
std::vector<double> vec_subtract(std::vector<double> const& x, std::vector<double> const& y)
{
	std::vector<double> diff (x.size());

	for(unsigned int i = 0; i < x.size(); i++)
	{
		diff[i] = x[i] - y[i];
	}
	return diff;
}

//do matrix addition between two vectors
std::vector<double> vec_add(std::vector<double> const& x, std::vector<double> const& y)
{
	std::vector<double> sum (x.size());

	for(unsigned int i = 0; i < x.size(); i++)
	{
		sum[i] = x[i] + y[i];
	}
	return sum;
}

//return vector product of CSR matrix and vector
std::vector<double> matvec_mult(std::vector<double> const& val,
								std::vector<int> const& row_ptr,
								std::vector<int> const& col_idx,
								std::vector<double> const& x)
{
	std::vector<double> prod (x.size());

	//index of product should match corresponding row
	//in csr, col index matches value index
	//go across each row, delineated by row_ptr
	for(unsigned int i = 0; i < x.size(); i++)
	{
		for(int j = row_ptr[i]; j < row_ptr[i+1]; j++)
		{
			prod[i] += val[j] * x[col_idx[j]];
		}
	}

	return prod;
}

//return vector product of CSR matrix and constant
void matconst_mult(std::vector<double>& val, std::vector<int> const& row_ptr,
					std::vector<int> const& col_idx, double const conval)
{
	for (unsigned int i = 0; i < val.size(); i++)
	{
		val[i] *= conval;
	}

}

//return product of constant and a vector
std::vector<double> constvec_mult(double const con, std::vector<double> const& x)
{
	std::vector<double> prod (x.size());

	for (unsigned int i = 0; i < x.size(); i++)
	{
		prod[i] = con * x[i];
	}

	return prod;
}

//return dot product of two vectors, transposing one
double dot_prod(std::vector<double> const& x, std::vector<double> const& y)
{
	double summing = 0.0;

	for (unsigned int i = 0; i < x.size(); i++)
	    {
			summing += x[i] * y[i];
	    }

	return summing;
}

//calculate the l2 norm of a vector
//vector const not changing
double L2norm(std::vector<double> const& x)
{
    double summing = 0.0;
    for (unsigned int i = 0; i < x.size(); i++)
    {
    	summing += x[i] * x[i];
    }
    return sqrt(summing);
}




