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


/*
   Pollard P-1 Factorization method.

   Algorithm 5.4.1 in "Prime Numbers, a Computational Perspective", 2ed. 
   R. Crandall, C. Pomerance.
*/

#include <vector>
#include <factorlab.hpp>

class base_item
{
   public:
   NTL::ZZ p;
   unsigned long a;
};

// {{{ FL::pollard_pm1()
NTL::ZZ FL::pollard_pm1(NTL::ZZ &n, NTL::ZZ &B)
{ 
   std::vector<base_item> base;
   NTL::ZZ k;
   NTL::ZZ c;
   NTL::ZZ g;
   NTL::ZZ p = NTL::to_ZZ(2);

   // Establish prime-power base
   while(p<B)
   {
      long a=0;
      k=p;
      while(k<B)
      {
         k*=p;
         a++;
      }

      base_item item;
      item.p = p;
      item.a = a;
      base.push_back(item);

      p = NTL::NextPrime(p+1);
   }

   // Perform power ladders
   c = 2;
   for(unsigned long i=0; i<base.size(); i++)
      for(unsigned long j=0; j<base[i].a; j++)
         c = NTL::PowerMod(c, base[i].p, n);
   
   // Test gcd
   g = NTL::GCD(c-1, n);

   return g;
}
// }}}


