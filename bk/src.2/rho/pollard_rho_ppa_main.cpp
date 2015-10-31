

#include <factorlab.hpp>

int main(int argc, char *argv[])
{
   if(argc < 4)
   {
      std::cout << "Usage: " << argv[0] << " <N> <k> <a>" << std::endl;
      exit(0);
   }

   NTL::ZZ n = NTL::to_ZZ(argv[1]);
   NTL::ZZ k = NTL::to_ZZ(argv[2]);
   NTL::ZZ a = NTL::to_ZZ(argv[3]);

   NTL::ZZ seed = NTL::to_ZZ(time(NULL));
   NTL::SetSeed(seed);

   std::cout << "Factoring " << n << " ..." << std::endl;

   if(NTL::ProbPrime(n, 10))
   {
      std::cout << "Prime: " << n << std::endl;
      return 0;
   }

   NTL::ZZ factor = FL::pollard_rho_ppa(n, k, a);

   std::cout << "Factor: " << factor << std::endl;

   return 0;
}



