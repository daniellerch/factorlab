
#include <iostream>
#include <fstream>

#include <factorlib.hpp>
#include <polynomial_selection.hpp>


// {{{ GNFS::get_base_m_expansion()
NTL::ZZX GNFS::get_base_m_expansion(int degree, NTL::ZZ &m, NTL::ZZ &n)
{
   NTL::ZZX f;
   NTL::ZZ N;
   
   f.rep.SetLength(degree);

   N = n-NTL::power(m, degree);
   for(int i=degree-1; i>0; i--) 
   {
      f.rep[i] = N / NTL::power(m, i);
      N -= NTL::power(m, i) * f.rep[i];;
   }


   f.rep[0] = N;
   for(int i=0; i<degree-1; i++) 
   {
      if(f.rep[i] > m/2) 
      {
         f.rep[i] -= m;
         f.rep[i+1] ++;
      }
   }

   return f;
}
// }}}

// {{{ GNFS::polynomial_selection()
void GNFS::polynomial_selection(GNFS::Polynomial &polynomial, Target &target)
{
   NTL::ZZ numZ;

   numZ = polynomial.d;
   polynomial.m = NTL::to_ZZ(NTL::pow(NTL::to_RR(target.n),1/NTL::to_RR(numZ)));
   polynomial.f = get_base_m_expansion(polynomial.d, polynomial.m, target.n);
  
   // f(x) reducible?
   if(FL::ZZX_is_reducible(polynomial.f, target.n))
   {
      std::cout << "\tf(x) is reducible: " << polynomial.f << std::endl;
      exit(1);
   }
   
}
// }}}





