#include <iostream>
#include <vector>

#include "COO2CSR.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"


	/* Method to modify sparse matrix dimensions */
void SparseMatrix::Resize(const int nrows, const int ncols)
{
	//the arguments that we pass in are divisions
	//add 1 to capture indices
	this->i_idx.resize(nrows);
	this->j_idx.resize(ncols);
	this->ncols = ncols;
	this->nrows = nrows;
}

/* Method to add entry to matrix in COO format */
void  SparseMatrix::AddEntry(const int i, const int j, const double val)
{
	a.push_back(val);
	i_idx.push_back(i);
	j_idx.push_back(j);
}

/* Method to convert COO matrix to CSR format using provided function */
void  SparseMatrix::ConvertToCSR()
{
	//call to method with current class instance
	COO2CSR(this->a, this->i_idx, this->j_idx);
}

/* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
std::vector<double> SparseMatrix::MulVec(std::vector<double> const& vec)
{
	std::cout << a.size() << " " << i_idx.size() << " " << j_idx.size() << " ";
	//first convert to CSR format
	this->ConvertToCSR();

	std::cout << a.size() << " " << i_idx.size() << " " << j_idx.size() << " ";

	//then use function, multiply matrix by arg vec, return
	return matvec_mult(this->a, this->i_idx, this->j_idx, vec);
}

/* Method to return the dimensions of the sparse matrix*/
std::vector<int> SparseMatrix::getDims() const
{
	std::vector<int> dims = { this->nrows, this->ncols };
	return dims;
}
/* TODO: Add any additional public methods you need */
