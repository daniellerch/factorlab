
/*
   Pollard P-1 Factorization method.

   Algorithm 5.4.1 in "Prime Numbers, a Computational Perspective", 2ed. 
   R. Crandall, C. Pomerance.
*/

#include <vector>
#include <factorlib.hpp>

class base_item
{
   public:
   NTL::ZZ p;
   unsigned long a;
};

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



