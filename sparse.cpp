#include <vector>

#include "COO2CSR.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"

class SparseMatrix
{
private:
	std::vector<int> i_idx;
	std::vector<int> j_idx;
	std::vector<double> a;

	std::vector<int>::iterator iter;

	int ncols;
	int nrows;

	/* TODO: Add any additional private data attributes and/or methods you need */


public:
	/* Method to modify sparse matrix dimensions */
	void Resize(int nrows, int ncols)
	{
		//the arguments that we pass in are divisions
		//add 1 to capture indices
		this->nrows = nrows;
		this->ncols = ncols;
		this->i_idx.resize(nrows+1);
		this->j_idx.resize(ncols+1);
	}

	/* Method to add entry to matrix in COO format */
	void  AddEntry(int i, int j, double val)
	{
		iter = i_idx.begin();
		i_idx.insert;
	}

	/* Method to convert COO matrix to CSR format using provided function */
	void  ConvertToCSR()
	{
		//call to method with current class instance
		COO2CSR(this->a, this->i_idx, this->j_idx);
	}

	/* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
	std::vector<double> MulVec(std::vector<double>& vec)
	{

		//first convert to CSR format
		this->ConvertToCSR();
		//then use function, multiply matrix by arg vec, return
		return matvec_mult(this->a, this->i_idx, this->j_idx, vec);
	}

	/* TODO: Add any additional public methods you need */

};