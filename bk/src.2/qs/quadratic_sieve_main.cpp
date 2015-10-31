

#include <factorlab.hpp>

int main(int argc, char *argv[])
{
   if(argc < 2)
   {
      std::cout << "Usage: " << argv[0] << " <N>" << std::endl;
      exit(0);
   }

   NTL::ZZ n = NTL::to_ZZ(argv[1]);

   NTL::ZZ seed = NTL::to_ZZ(time(NULL));
   NTL::SetSeed(seed);

   std::cout << "Factoring " << n << " ..." << std::endl;

   NTL::ZZ factor = FL::quadratic_sieve(n);
   std::cout << "Factor: " << factor << std::endl;

   return 0;
}



