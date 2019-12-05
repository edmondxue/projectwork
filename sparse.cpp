#include <iostream>
#include <vector>

#include "COO2CSR.hpp"
#include "matvecops.hpp"
#include "sparse.hpp"


	/* Method to modify sparse matrix dimensions */
void SparseMatrix::Resize(const int nrows, const int ncols)
{
	this->ncols = ncols;
	this->nrows = nrows;
}

/* Method to add entry to matrix in COO format */
void  SparseMatrix::AddEntry(const int i, const int j, const double val)
{
	a.push_back(val);
	i_idx.push_back(i);
	j_idx.push_back(j);

	//if (j_idx.size() > 7200)
	//{ 
	//	std::cout << j_idx.size() << " ";
	//	std::cout << i_idx.size() << " \n";
	//}
		
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
	//std::cout << "BEFORE" << "\n ";
	//std::cout << a.size() << " \n";
	//std::cout << i_idx.size() << " \n";
	//std::cout << j_idx.size() << " \n";

	//for (auto i = a.begin(); i != a.end(); ++i)
	//	std::cout << *i << ' ';

	//std::cout << "\n";
	//for (auto i = i_idx.begin(); i != i_idx.end(); ++i)
	//	std::cout << *i << ' ';
	//std::cout << "\n";
	//for (auto i = j_idx.begin(); i != j_idx.end(); ++i)
	//	std::cout << *i << ' ';


	//std::cout << a.size() << " \n";
	//std::cout << i_idx.size() << " \n";
	//std::cout << j_idx.size() << " \n";

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
