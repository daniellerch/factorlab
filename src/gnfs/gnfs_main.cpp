
#include <gnfs.hpp>
#include <polynomial_selection.hpp>
#include <factor_base.hpp>
#include <sieve.hpp>
#include <linear_algebra.hpp>
#include <square_root.hpp>


// {{{ calc_B()
// Prime Numbers: A Computational Perspective - Crandall & Pomerance
// B = exp( (8/9)^(1/3) * (ln n)^(1/3) * (lnln n)^(2/3) )
int
calc_B(NTL::ZZ n)
{
   NTL::RR lnn;
   NTL::RR lnlnn;
   NTL::RR A, B, C;
   int result = 0;

   lnn = NTL::log(n);
   lnlnn = NTL::log(lnn);

   A = NTL::pow(NTL::to_RR((double)8 / 9), NTL::to_RR((double)1 / 3));
   B = NTL::pow(lnn, NTL::to_RR((double)1 / 3));
   C = NTL::pow(lnlnn, NTL::to_RR((double)2 / 3));

   result = 10*NTL::to_int(exp(A * B * C));

   return result;
}

// }}}

// {{{ calc_U()
// Factoring Integers With The Number Field Sieve 
// Buhler & Lenstra & Pomerance
// u = exp( 1/2 * ( d*(lnd)+sqrt( (d*lnd)^2 + 4*(ln(n^(1/d)))*lnln(n^(1/d))) ))
int
calc_U(NTL::ZZ & n, int d)
{
   int result = 0;

   NTL::RR lnd;
   NTL::RR lnn1d;
   NTL::RR lnlnn1d;
   NTL::RR A, B, C, D;

   lnd = NTL::log(NTL::to_RR(d));
   lnn1d = NTL::log(NTL::pow(NTL::to_RR(n), NTL::to_RR((double)1 / d)));
   lnlnn1d = NTL::log(lnn1d);

   double e = 0;

   A = (double)((1 + e) / 2);
   B = d * lnd;
   C = (d * lnd) * (d * lnd);
   D = 4 * lnd * lnlnn1d;

   result = to_int(NTL::exp(A * (B + NTL::sqrt(C + D))));

   return result;
}

// }}}

// ----------------------------------------------------------------------------
// MAIN - General Number Field Sieve
// ----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   NTL::ZZ seed = NTL::to_ZZ(time(NULL));
   NTL::SetSeed(seed);

   int t0, t1;
   int u;
   int v;
   int pairs_needed;
   NTL::ZZ xZ;
   NTL::ZZ yZ;      
   std::ifstream fin;
   std::string line = "  -------------------------------------"
                      "-------------------------------------  ";

	GNFS::Target target;


   if(argc != 2)
   {
      std::cout << "Usage: " << argv[0] << " [N] " << std::endl << std::endl;;
      exit(1);
   }

   target.n = NTL::to_ZZ(argv[1]);
	target.nbits = NTL::NumBits(target.n);
   NTL::RR::SetPrecision(NTL::NumBits(target.n));

	// Create polynomial
	GNFS::Polynomial polynomial;
	polynomial.d=3;


   target.digits = (int)(log(target.n) / log(10)) + 1;
   double div=1;
   //int div=1;
   target.t = calc_U(target.n, polynomial.d)/div;
   target.C = calc_B(target.n)/div;
   
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tGeneral Number Field Sieve (GNFS)" << std::endl;
   std::cout << line << std::endl;

   std::cout << "\tTarget Number: " << target.n << std::endl;
   std::cout << "\tDigits: " << target.digits << std::endl;
   std::cout << "\tNum Bits: " << target.nbits << std::endl;
   std::cout << "\tdegree:  " << polynomial.d << std::endl;
   std::cout << "\tFB Size: " << target.t << std::endl;
   std::cout << "\tSive Interval: " << target.C << std::endl;


   // ------------------------------------------------------------------------
   // STEP 1: Polynomial Selection
   // ------------------------------------------------------------------------
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tPolynomial Selection" << std::endl;
   std::cout << line << std::endl;
   
   t0 = time(NULL);

   polynomial_selection(polynomial, target);

	std::cout << "\tPolynomial: " << polynomial.f << std::endl;
	std::cout << "\tDegree: " << polynomial.d << std::endl;
	std::cout << "\tm: " << polynomial.m << std::endl;

   t1 = time(NULL);
   printf("\ttime: %d\n", t1-t0);



   // ------------------------------------------------------------------------
   // STEP 2: Create Factor Base
   // ------------------------------------------------------------------------
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tMake Factor Base" << std::endl;
   std::cout << line << std::endl;

   t0 = time(NULL);

   FactorBase fb;

   fb.make_RFB(polynomial, target);
   std::cout << "\tRFB: " << fb.RFB.size() << " elements" << std::endl;

   fb.make_AFB(polynomial, target);
   std::cout << "\tAFB: " << fb.AFB.size() << " elements" << std::endl;

   fb.make_QFB(target, polynomial, fb.AFB[fb.AFB.size()-1]);
   std::cout << "\tQCB: " << fb.QCB.size() << " elements" << std::endl;

   v = target.digits;
   u = polynomial.d * target.t;

   t1 = time(NULL);
   printf("\ttime: %d\n", t1-t0);


   // ------------------------------------------------------------------------
   // STEP 3: Sieve
   // ------------------------------------------------------------------------
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tSieve" << std::endl;
   std::cout << line << std::endl;

   t0 = time(NULL);

   std::vector<int> av; 
   std::vector<int> bv; 
   pairs_needed = target.t + u + v + 2;
   sieve(polynomial, target, fb, pairs_needed, av, bv);

   t1 = time(NULL);
   printf("\n\n\ttime: %d\n", t1-t0);
 

   // ------------------------------------------------------------------------
   // STEP 4: Linear Algebra
   // ------------------------------------------------------------------------
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tLinear Algebra" << std::endl;
   std::cout << line << std::endl;
   t0 = time(NULL);

   Matrix matrix(pairs_needed, pairs_needed);
   linear_algebra(polynomial, target, fb, matrix, av, bv);

   t1 = time(NULL);
   printf("\n\n\ttime: %d\n", t1-t0);

	
   // ------------------------------------------------------------------------
   // STEP 5: Square Root
   // ------------------------------------------------------------------------
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tSquare Root" << std::endl;
   std::cout << line << std::endl;

   t0 = time(NULL);

   square_root(polynomial, target, matrix, pairs_needed, pairs_needed-1, 
      fb, av, bv, xZ, yZ);

   t1 = time(NULL);
   printf("\ttime: %d\n", t1-t0);


   // ------------------------------------------------------------------------
   // STEP 6: Find Factors
   // ------------------------------------------------------------------------
   std::cout << std::endl;
   std::cout << line << std::endl;
   std::cout << "\tFactors: " << std::endl;
   std::cout << line << std::endl;

   std::cout << "\tfactor: " << NTL::GCD(xZ-yZ, target.n) << std::endl;
   std::cout << "\tfactor: " << NTL::GCD(xZ+yZ, target.n) << std::endl;
   std::cout << std::endl;

   return (EXIT_SUCCESS);
}


