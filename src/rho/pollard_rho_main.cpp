

#include <factorlib.hpp>

int main(int argc, char *argv[])
{
   if(argc < 2)
   {
      std::cout << "Usage: " << argv[0] << " <N> <k(optional)>" << std::endl;
      exit(0);
   }

   NTL::ZZ n = NTL::to_ZZ(argv[1]);

   long k = 2;
   if(argc==3)
   {
      k=atol(argv[2]);
      std::cout << "Using k=" << k << std::endl;
   }

   NTL::ZZ seed = NTL::to_ZZ(time(NULL));
   NTL::SetSeed(seed);

   std::cout << "Factoring " << n << " ..." << std::endl;

   if(NTL::ProbPrime(n, 10))
   {
      std::cout << "Prime: " << n << std::endl;
      return 0;
   }

   NTL::ZZ factor = FL::pollard_rho(n, k);

   std::cout << "Factor: " << factor << std::endl;

   return 0;
}



