
#include <iostream>
#include <fstream>
#include <polynomial_selection.hpp>
#include <factorlab.hpp>


// {{{ FL::ZZX_is_reducible()
bool FL::ZZX_is_reducible(NTL::ZZX &f, NTL::ZZ &n)
{
   NTL::ZZ div, root, val, product;

   root = NTL::SqrRoot(abs(f.rep[0])) + 1;
   for(div = 1; div < root; div++)
   {
      if(n % div == 0)
      {
         val = div;
         product = 1;
         product *= FL::ZZX_evaluate(f, val);
         val = 0 - div;
         product *= FL::ZZX_evaluate(f, val);
         val = n / div;
         product *= FL::ZZX_evaluate(f, val);
         val = 0 - (n / div);
         product *= FL::ZZX_evaluate(f, val);
         if(product == 0)
         {
            return true;
         }
      }
   }
   return false;
}
// }}}

// {{{ FL::ZZX_evaluate()
// Polynomial function
NTL::ZZ FL::ZZX_evaluate(const NTL::ZZX &f, NTL::ZZ &x)
{
   NTL::ZZ temp;

   temp = power(x, NTL::deg(f)+1);
   for(int i = 0; i < NTL::deg(f)+1; i++)
      temp += NTL::coeff(f, i) * NTL::power(x, i);

   return temp;
}
// }}}

// {{{ FL::ZZX_evaluate()
NTL::RR FL::ZZX_evaluate(const NTL::ZZX &f, NTL::RR &x)
{
   NTL::RR temp;

   temp = power(x, NTL::deg(f)+1);
   for(int i = 0; i < NTL::deg(f)+1; i++)
      temp += NTL::to_RR(NTL::coeff(f, i)) * NTL::to_RR(NTL::power(x, i));

   return temp;
}
// }}}





