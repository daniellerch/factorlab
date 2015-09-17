
#ifndef __GNFS_HPP__
#define __GNFS_HPP__

#include <iostream>
#include <fstream>
#include <vector>

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_GF2.h>
#include <NTL/ZZX.h>

#include <NTL/GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>

namespace GNFS
{
   class Target
   {
      public:
         NTL::ZZ n;
         int digits;
         int nbits;
         int t;
         int C;
   };
}

#endif



