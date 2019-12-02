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

	/* TODO: Add any additional private data attributes and/or methods you need */


public:
	/* Method to modify sparse matrix dimensions */
	void Resize(const int nrows, const int ncols)
	{
		//the arguments that we pass in are divisions
		//add 1 to capture indices
		this->i_idx.resize(nrows+1);
		this->j_idx.resize(ncols+1);
	}

	/* Method to add entry to matrix in COO format */
	void  AddEntry(const int i, const int j, const double val)
	{
		a.push_back(val);
		i_idx.push_back(i);
		j_idx.push_back(j);
	}

	/* Method to convert COO matrix to CSR format using provided function */
	void  ConvertToCSR()
	{
		//call to method with current class instance
		COO2CSR(this->a, this->i_idx, this->j_idx);
	}

	/* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
	std::vector<double> MulVec(const std::vector<double> const& vec)
	{
		//first convert to CSR format
		this->ConvertToCSR();
		//then use function, multiply matrix by arg vec, return
		return matvec_mult(this->a, this->i_idx, this->j_idx, vec);
	}

	/* TODO: Add any additional public methods you need */

};