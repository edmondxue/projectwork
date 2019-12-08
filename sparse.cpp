#include <vector>

#include "COO2CSR.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"


	/* Method to modify sparse matrix dimensions */
void SparseMatrix::Resize(const int nrows, const int ncols)
{
	this->nrows = nrows;
	this->ncols = ncols;
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
	//use function, multiply matrix by arg vec, return

	std::vector <double> prod;
	prod = matvec_mult(this->a, this->i_idx, this->j_idx, vec);

	return prod;
}

/* Method to perform sparse matrix constant multiplication using CSR formatted matrix */
void SparseMatrix::MulConst(double const conval)
{

	//then use function, multiply matrix by arg vec, return
	matconst_mult(this->a, this->i_idx, this->j_idx, conval);
}

/* Method to return the dimensions of the sparse matrix*/
std::vector<int> SparseMatrix::getDims() const
{
	std::vector<int> dims = { this->nrows, this->ncols };
	return dims;
}
