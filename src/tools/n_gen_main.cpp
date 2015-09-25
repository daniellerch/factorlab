
#include <iostream>
#include <NTL/ZZ.h>

int main(int argc, char *argv[])
{
   if(argc!=2)
   {
      std::cout << "Usage: " << argv[0] << " [Num Bits] " << std::endl;
      return 0;
   }

   NTL::ZZ p, q, n;
   long numbits = atol(argv[1]);

   NTL::SetSeed(NTL::to_ZZ(time(NULL)));
   p = NTL::RandomPrime_ZZ(numbits);
   q = NTL::RandomPrime_ZZ(numbits);
   n = p*q;

   std::cout << "p: " << p << " * " << std::endl;    
   std::cout << "q: " << q << " = " << std::endl;    
   std::cout << "n: " << n << std::endl;    

   return 0;
}




