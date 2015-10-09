/*
   Copyright 2015 Daniel Lerch Hostalot.

   This file is part of factorlab.

   factorlab is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   FLINT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FLINT; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */



#include <ecm/ecm.hpp>

int main(int argc, char *argv[])
{
   if(argc < 2)
   {
      std::cout << "Usage: " << argv[0] << " <N> [B1] [B2]" << std::endl;
      exit(0);
   }

   NTL::ZZ seed = NTL::to_ZZ(time(NULL));
   NTL::SetSeed(seed);

   NTL::ZZ n = NTL::to_ZZ(argv[1]);
   long B1=0;
   long B2=0;

   if(argc>=3)
      B1=atol(argv[2]);

   if(argc>=4)
      B2=atol(argv[3]);

   std::cout << "Factoring " << n << " ..." << std::endl;

   NTL::ZZ factor = ecm(n, B1, B2);
   std::cout << "Factor: " << factor << std::endl;

   return 0;
}



