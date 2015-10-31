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



#include <vector>
#include <factorlab.hpp>

NTL::ZZ F1(NTL::ZZ &x, NTL::ZZ &n)
{
   // f(x) = ((x^2 - 85)^2 - 4176)^2 - 2880^2

   NTL::ZZ fx;
   fx = NTL::PowerMod(NTL::PowerMod(x, 2, n) - 85, 2, n);
   fx = NTL::PowerMod(fx - 4176, 2, n);
   fx = fx - 8294400;

   return fx;
}

NTL::ZZ F2(NTL::ZZ &x, NTL::ZZ &n)
{
   // f(x) = (((x^2 - 67405)^2 - 3525798096)^2 -  533470702551552000)^2 - 469208209191321600^2

   NTL::ZZ fx;
   fx = NTL::PowerMod(NTL::PowerMod(x, 2, n) - 67405, 2, n);
   fx = NTL::PowerMod(fx - 3525798096, 2, n);
   fx = NTL::PowerMod(fx - 533470702551552000, 2, n);
   fx = fx - NTL::to_ZZ("220156343572527011594632754626560000");

   return fx;
}


NTL::ZZ F(NTL::ZZ &x, NTL::ZZ &n, long begin, long block)
{
   NTL::ZZ r = NTL::to_ZZ(1);

   for(long i=begin; i<begin+block; i++)
      r *= x - i;

   return r;
}

// {{{ FL::polynomial_evaluacion()
NTL::ZZ FL::polynomial_evaluation(NTL::ZZ &n)
{
   NTL::ZZ g; 
   NTL::ZZ x; 

   long block=100;
   unsigned long cnt=0;
   long i=0;
   for(;;)
   {
      cnt++;
      x = NTL::RandomBnd(n);
      g = NTL::GCD( F(x, n, i, block), n);
      i+=block;
   
      if(g>1 && g<n)
         break;
   }
   
   std::cout << "Iterations: " << cnt << std::endl;

   return g;
}
// }}}


