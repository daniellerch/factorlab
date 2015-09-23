
#ifndef __GNFS_FACTOR_BASE_HPP__
#define __GNFS_FACTOR_BASE_HPP__


#include <iostream>
#include <fstream>

#include <factorlib.hpp>
#include <gnfs.hpp>
#include <polynomial_selection.hpp>


class FactorBase
{
   public:
   std::vector <int> RFB;
   std::vector <int> RFBm;
   std::vector <double> RFBlog;
   std::vector <int> AFB;
   std::vector <int> AFBr;
   std::vector <double> AFBlog;
   std::vector <int> QCB;
   std::vector <int> QCBs;


// {{{ prime_logarithm()
double prime_logarithm(std::vector <int> &FB, int bound, int p, NTL::ZZ n)
{
	NTL::ZZ x;
	double res=0;
	
	x = n-p;
	for(int i=0; i<bound; i++)
	{
		double l = FB[i];
		if(x%FB[i]==0)
			res += log(l);
		else
			return 0;
	}

	return res;
}
// }}}

   // {{{ make_RFB()
   void make_RFB(GNFS::Polynomial &polynomial, GNFS::Target &target)
   {
		int bound = target.t;
		int ext = 0;
      long p=0;

      for(int i = 0; i < bound; i++)
      {
         p = NTL::NextPrime(p+1);
			double l;

			l = p;
         RFB.push_back(p);
         RFBm.push_back(polynomial.m % p);

			if(i < bound-ext)
				l = log(l);
			else
				l = prime_logarithm(RFB, bound-ext, p, target.n);

         RFBlog.push_back(l);
      }
   }
   // }}}

    // {{{ make_AFB()
   void make_AFB(GNFS::Polynomial &polynomial, GNFS::Target &target)
   {
      int u = 0;
      long p=0;
		double l;
		int bound = target.t;
		int ext = 0;

      NTL::ZZ jZ, pZ;

      while(u < polynomial.d * bound)
      {
         p = NTL::NextPrime(p+1);
         pZ = p;
         int count = 0;

         for(int j = 0; j < p && count < polynomial.d; j++)
         {
            jZ = j;
            if(FL::ZZX_evaluate(polynomial.f, jZ) % pZ == 0)
            {
               count++;
               AFB.push_back(p);
               AFBr.push_back(j);

					if(u < polynomial.d*bound-ext)
						l = log(p);
					else
						l = prime_logarithm(AFB, bound-ext, p, target.n);

               AFBlog.push_back(l);
               u++;
            }
         }
      }
   }

   // }}}

   // {{{ make_QFB()
   void make_QFB(GNFS::Target &target, GNFS::Polynomial &polynomial, int lastp)
   {
      int p = lastp;
      int v = 0;
      NTL::ZZ jZ, pZ;

      while(v<target.digits)
      {
         p = NTL::NextPrime(p+1);

         pZ = p;
         for(int j = 0; j < p; j++)
         {
            jZ = j;

            if(FL::ZZX_evaluate(polynomial.f, jZ) % pZ == 0)
            {
               QCB.push_back(p);
               QCBs.push_back(j);
               v++;
            }
         }
      }
   }

   // }}}

};

#endif

