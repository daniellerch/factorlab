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

// {{{ create_factor_base()
static void create_factor_base(NTL::ZZ &n)
{


}
// }}}

// {{{ FL::quadratic_sieve()
NTL::ZZ FL::quadratic_sieve(NTL::ZZ &n)
{
   NTL::ZZ B = NTL::to_ZZ(10000000);

   std::cout << "B: " << B << std::endl;


   create_factor_base(n);

   // TODO

   std::cout << "Warning: not implemented!" << std::endl;
   return n;   
}
// }}}


