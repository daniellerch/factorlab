
#ifndef __GGNFS_NTL_EXTENSION_HPP__
#define __GGNFS_NTL_EXTENSION_HPP__

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>

namespace NTL_extension
{

   bool ZZX_is_reducible(NTL::ZZX &f, NTL::ZZ &n);

   NTL::ZZ ZZX_evaluate(const NTL::ZZX &f, NTL::ZZ &x);

   NTL::RR ZZX_evaluate(const NTL::ZZX &f, NTL::RR &x);

}

#endif

