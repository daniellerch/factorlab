
#ifndef __GGNFS_POLYNOMIAL_SELECTION_HPP__
#define __GGNFS_POLYNOMIAL_SELECTION_HPP__


#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>

#include <gnfs.hpp>

namespace GNFS
{

   class Polynomial
   {
      public:
         NTL::ZZX f;
         NTL::ZZ m;
         int d;

   };


   NTL::ZZ dF(NTL::ZZX &poly, NTL::ZZ &x);

   NTL::ZZX get_base_m_expansion(int degree, NTL::ZZ &m, NTL::ZZ &n);

   void polynomial_selection(GNFS::Polynomial &polynomial, Target &target);

}

#endif

