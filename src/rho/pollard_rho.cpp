/*
   Copyright 2015 Daniel Lerch Hostalot.

   This file is part of FactorLib.

   FactorLib is free software; you can redistribute it and/or modify                
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
   Pollard Rho Factorization method.

   Algorithm 5.2.1 in "Prime Numbers, a Computational Perspective", 2ed. 
   R. Crandall, C. Pomerance.
*/

#include <factorlib.hpp>

// {{{ F()
static inline NTL::ZZ F(NTL::ZZ &x, NTL::ZZ &a, long k, NTL::ZZ &n)
{
   // f(x) = x^2k + a (mod n)
   return NTL::AddMod(NTL::PowerMod(x, k, n), a, n);
}
// }}}

// {{{ FL::pollard_rho()
NTL::ZZ FL::pollard_rho(NTL::ZZ &n, long k)
{   
   // Usually we use k=2 with complexity c*sqrt(p) but it is an established 
   // heuristic that the expected number of iterations can be reduced from 
   // c*sqrt(p) to c*sqrt(p)/sqrt(gcd(p-1,2k)-1). 

   NTL::ZZ a;
   NTL::ZZ s;
   NTL::ZZ U;
   NTL::ZZ V;
   NTL::ZZ g;

   unsigned long cnt=0;
   for(;;)
   {
      // Choose seeds
      a=NTL::RandomBnd(n-3)+1;
      s=NTL::RandomBnd(n);
      std::cout << "Seed: " << a << ", " << s << std::endl;
      U=s;
      V=s;

      /*
      int rnd_points=1000;
      NTL::ZZ Vr[rnd_points];
      for(int i=0; i<rnd_points; i++)
         Vr[i]=NTL::RandomBnd(n);
      */

      //bool found=false;
      for(;;)
      {
         cnt++;

         // Factor search
         U = F(U, a, k, n);
         V = F(V, a, k, n);
         V = F(V, a, k, n);
 

         /*
         for(int i=0; i<rnd_points && !found; i++)
         {
            Vr[i] = F(Vr[i], a, k, n);
            Vr[i] = F(Vr[i], a, k, n);
            g = NTL::GCD(U-Vr[i], n);
            if(1<g && g<n)
            {
               std::cout << "Factor(" << i << "): " << g << ", it: " << cnt << std::endl;
               found = true;
            }
         }
         */
   
         g = NTL::GCD(U-V, n);

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


