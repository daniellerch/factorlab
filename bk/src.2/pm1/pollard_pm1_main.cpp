

#include <factorlab.hpp>

int main(int argc, char *argv[])
{
   if(argc < 3)
   {
      std::cout << "Usage: " << argv[0] << " <N> <B>" << std::endl;
      exit(0);
   }

   NTL::ZZ n = NTL::to_ZZ(argv[1]);
   NTL::ZZ B = NTL::to_ZZ(argv[2]);


   NTL::ZZ seed = NTL::to_ZZ(time(NULL));
   NTL::SetSeed(seed);

   std::cout << "Factoring " << n << " ..." << std::endl;

   NTL::ZZ factor = FL::pollard_pm1(n, B);
   std::cout << "Factor: " << factor << std::endl;

   return 0;
}



