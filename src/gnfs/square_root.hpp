
#ifndef __GNFS_SQUARE_ROOT_HPP__
#define __GNFS_SQUARE_ROOT_HPP__

#include <gnfs.hpp>
#include <factor_base.hpp>
#include <sieve.hpp>
#include <linear_algebra.hpp>

void square_root(GNFS::Polynomial &polynomial, GNFS::Target &target, 
	Matrix &matrix, int pairs_needed, int nRows,  
	FactorBase &fb,  const std::vector<int> &av, const std::vector<int> &bv,
   NTL::ZZ &xZ, NTL::ZZ &yZ);

void find_roots(GNFS::Polynomial &polynomial);

NTL::RR Newton(GNFS::Polynomial &polynomial, NTL::RR start);

bool has_roots(GNFS::Polynomial &polynomial);



#endif

