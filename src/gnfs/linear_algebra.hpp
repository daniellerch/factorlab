
#ifndef __GNFS_LINEAR_ALGEBRA_HPP__
#define __GNFS_LINEAR_ALGEBRA_HPP__

#include <gnfs.hpp>
#include <factor_base.hpp>
#include <matrix/smat_long.h>

class Matrix
{
   public:
   NTL::smat_long sM;
   //NTL::mat_GF2 sM;
   NTL::svec_long sdependent;
   //NTL::vec_GF2 sfreeCols;
   NTL::svec_long sfreeCols;

   Matrix(int col, int row)
   {
		sM.SetDims(col, row);  
      sdependent.SetLength(row);
      sfreeCols.SetLength(row);
   }

};


void linear_algebra(GNFS::Polynomial &polynomial, GNFS::Target &target, 
	FactorBase &fb, Matrix &matrix, const std::vector<int> &av, 
   const std::vector<int> &bv);

#endif

