
#ifndef __GNFS_LINEAR_ALGEBRA_HPP__
#define __GNFS_LINEAR_ALGEBRA_HPP__

#include <gnfs.hpp>
#include <factor_base.hpp>

class Matrix
{
   public:
   NTL::mat_GF2 sM;
   NTL::vec_GF2 sdependent;
   NTL::vec_GF2 sfreeCols;
   int col;
   int row;

   Matrix(int col, int row)
   {
      this->col = col;
      this->row = row;
		sM.SetDims(col, row);  
      sdependent.SetLength(row);
      sfreeCols.SetLength(row);
   }

};


void linear_algebra(GNFS::Polynomial &polynomial, GNFS::Target &target, 
	FactorBase &fb, Matrix &matrix, const std::vector<int> &av, 
   const std::vector<int> &bv);

#endif

