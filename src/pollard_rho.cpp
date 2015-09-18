
/*
   Pollard Rho Factorization method.

   Algorithm 5.2.1 in Prime Numbers, a Computational Perspective, 2ed. 
   R. Crandall, C. Pomerance.
*/

#include <factorlib.hpp>

inline NTL::ZZ F(NTL::ZZ &x, NTL::ZZ &a, long k, NTL::ZZ &n)
{
   return NTL::AddMod(NTL::PowerMod(x, k, n), a, n);
}

NTL::ZZ pollard_rho(NTL::ZZ &n, long k)
{   
   // Usually we use k=2 with complexity c*sqrt(p) but it is an established 
   // heuristic that the expected number of iterations can be reduced from 
   // c*sqrt(p) to c*sqrt(p)/sqrt(gcd(p-1,2k)-1). 

   NTL::ZZ a;
   NTL::ZZ s;
   NTL::ZZ U;
   NTL::ZZ V;
   NTL::ZZ g;

   for(;;)
   {
      // Choose seeds
      a=NTL::RandomBnd(n-3)+1;
      s=NTL::RandomBnd(n);
      
      U=s;
      V=s;

      for(;;)
      {
         // Factor search
         U = F(U, a, k, n);
         V = F(V, a, k, n);
         V = F(V, a, k, n);
         g = NTL::GCD(U-V, n);

         if(g==1) 
            continue;
         else if(g==n)
            break;

         // Success
         return g;
      }
   }
}



