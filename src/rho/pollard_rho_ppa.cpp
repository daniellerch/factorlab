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



#include <factorlab.hpp>

// {{{ F()
static inline NTL::ZZ F(NTL::ZZ &x, NTL::ZZ &k, NTL::ZZ &a, NTL::ZZ &b, NTL::ZZ &n)
{
   // x^p = x mod p
   // x^a * x^p = x^a *x mod p
   // x^a * x^p = x^a *x mod p
   // x^(p+a) = x^(a+1) mod p

   // k = p+a
   // f(x) = x^k - x^(a+1) + b (mod n)

   // k = m(p+a)
   // f(x) = x^mk - x^m(a+1) + b (mod n)
   
   // k == m(p+a) ~ B-smooth
   // f(x) = x^k - x^Bi(a+1) + b (mod n) for every subset of B. We can discard
   // a subset if it does not find factors in one iteration. (ufff!)


   NTL::ZZ c1 = NTL::PowerMod(x, k, n);
   NTL::ZZ c2 = NTL::PowerMod(x, a+1, n);
   NTL::ZZ c = NTL::AddMod(c1, -c2+n, n);

   return NTL::AddMod(c, b, n);
}
// }}}

// {{{ FL::pollard_rho_ppa()
NTL::ZZ FL::pollard_rho_ppa(NTL::ZZ &n, NTL::ZZ &k, NTL::ZZ &a)
{   
   // Usually we use k=2 with complexity c*sqrt(p) but it is an established 
   // heuristic that the expected number of iterations can be reduced from 
   // c*sqrt(p) to c*sqrt(p)/sqrt(gcd(p-1,2k)-1). 

   NTL::ZZ b;
   NTL::ZZ s;
   NTL::ZZ U;
   NTL::ZZ V;
   NTL::ZZ g;

   unsigned long cnt=0;
   for(;;)
   {
      // Choose seeds
      b=NTL::RandomBnd(n-3)+1;
      s=NTL::RandomBnd(n);
      //std::cout << "Seed: " << b << ", " << s << std::endl;
      U=s;
      V=s;

      for(;;)
      {
         cnt++;

         // Factor search
         U = F(U, k, a, b, n);
         V = F(V, k, a, b, n);
         V = F(V, k, a, b, n);
         g = NTL::GCD(U-V, n);
         
         //std::cout << "U: " << U << std::endl;
         //std::cout << "V: " << V << std::endl;
         //std::cout << "GCD: " << g << std::endl;
         
         if(g==1) 
            continue;
         else if(g==n)
            break;

         // Success
         std::cout << "Iterations: " << cnt << std::endl;
         return g;
      }
   }
}
// }}}


