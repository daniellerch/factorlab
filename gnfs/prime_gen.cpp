
#include <iostream>
#include <NTL/ZZ.h>

int main(int argc, char *argv[])
{
   if(argc!=2)
   {
      std::cout << "Usage: " << argv[0] << " [Max Prime] " << std::endl;
      return 0;
   }

   NTL::ZZ max = NTL::to_ZZ(argv[1]);
   NTL::ZZ prime = NTL::to_ZZ(2);

   std::cout << 2 << std::endl;
   while(prime<max)
   {
      prime = NextPrime(prime+1);
      std::cout << prime << std::endl;
   }
   return 0;
}




