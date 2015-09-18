
#include <sieve.hpp>

#define SKIP -321.123


// {{{ algebraic_norm()
//  NORM(a,b) = (-b)^d * f(-a/b)
NTL::ZZ algebraic_norm(GNFS::Polynomial &polynomial, int a, int b)
{
	// TODO: degree 3 especific
   NTL::ZZ temp, negOne, A, B;

   negOne = -1;
   A = a;
   B = b;
   temp = NTL::power(A, polynomial.d);
   for(int i = 0; i < polynomial.d; i++)
   {
      temp += NTL::power(negOne, i + 3) * 
			(polynomial.f.rep[i]) * NTL::power(A, i) * NTL::power(B, 3 - i);
   }
   return temp;
}
// }}}

// {{{ sieve()
void sieve(GNFS::Polynomial &polynomial, GNFS::Target &target, FactorBase &fb,   
   int pairs_needed, std::vector<int> &av, std::vector<int> &bv)
{
   std::vector<double> normv; 
   int b = 0;
   int pairs_found = 0;
   normv.resize(2*target.C+1);
   av.resize(pairs_needed);
   bv.resize(pairs_needed);
   int u = polynomial.d*target.t;
   int i;
   int p;
   NTL::ZZ numZ;
   NTL::ZZ valZ;


   while(pairs_found < pairs_needed)
   {
      // Note: i-C = a
      int a;

      b++;
			
      for(i=0; i<=2*target.C; i++)
      {
			valZ = (i-target.C) + b*polynomial.m;

         if(NTL::GCD(i-target.C, b)!=1 || valZ==0)
            normv[i] = SKIP;

         else
				normv[i] = 0.5 - log(abs(valZ));
      }


      // Test for RFB-smoothness
      for(i=0; i < target.t; i++)
      {
         p = fb.RFB[i];
         a = (int)(p * (-floor(target.C / p)) - (b * fb.RFBm[i])) + target.C;
         
			while(a < 0)
            a += p;
         
         while(a <= 2*target.C)
         {
            if(normv[a] != SKIP)
               normv[a] += fb.RFBlog[i];

            a += p;
         }

         // Now check p^2, p^3, ....
         int pPow = p * p;

         while(pPow < 2 * target.C)
         {
            a = (int)(pPow * (-floor(target.C/pPow))-
					(b*(polynomial.m % pPow))) + target.C;

            while(a < 0)
               a += pPow;
            
            while(a <= 2*target.C)
            {
               if(normv[a] != SKIP)
                  normv[a] += fb.RFBlog[i];
               
               a += pPow;
            }
            pPow *= p;
         }
      }

      // Initialize normv with -ln( (-b)^d*f(-a/b) )
      for(i=0; i<=2*target.C; i++)
      {
         if(normv[i] > 0)
         {
            // init Norm
            numZ = abs(algebraic_norm(polynomial, i-target.C, b));
            if(numZ == 0)
               normv[i] = SKIP;
            else
               normv[i] = fb.AFBlog[u-1] - log(numZ);
         }
         else
            normv[i] = SKIP;
      }

      // Test for AFB-smoothness
      for(i=0; i<u; i++)
      {
         p = fb.AFB[i];

         a = (int)(p * (-floor(target.C / p)) - (b * fb.AFBr[i])) + target.C;
         while(a < 0)
            a += p;

         while(a <= 2*target.C)
         {
            if(normv[a] != SKIP)
               normv[a] += fb.AFBlog[i];

            a += p;
         }
      }

      // Search for good (a,b)'s
      for(int i=0; i<=2*target.C; i++)
      {
         if(normv[i] > 0)
         {
            if(pairs_found < pairs_needed)
            {
               av[pairs_found] = i-target.C;
               bv[pairs_found++] = b;
            }
         }
      }

		// info
      std::cout << "\tit: " << b << " relation: " << pairs_found
			<< " needed: " << pairs_needed << " completed: (" 
			<< (pairs_found*100)/pairs_needed << "%)\r";
   }


}
// }}}


