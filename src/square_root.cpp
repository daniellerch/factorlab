
#include <NTL/ZZX.h>
#include <polynomial_selection.hpp>
#include <square_root.hpp>


// {{{ Newton()
NTL::RR Newton(GNFS::Polynomial &polynomial, NTL::RR start)
{
   NTL::RR error = NTL::to_RR(0.001);
	NTL::RR x0;
	NTL::RR y0;
	NTL::RR yp;
	NTL::RR x1;

   x0 = 0;
   x1 = start;

	while(NTL_extension::ZZX_evaluate(polynomial.f, x1) > error)
	{
      x0 = x1;

      y0=NTL_extension::ZZX_evaluate(polynomial.f, x0);
      yp=NTL_extension::ZZX_evaluate(NTL::diff(polynomial.f), x0);
      x1=x0-y0/yp;
   }

   return x1;
}
// }}}

// {{{ dF()
NTL::ZZ dF(NTL::ZZX &poly, NTL::ZZ &x)
{
   // TODO: degree 3 specific
   NTL::ZZ temp;

   temp = 3 * NTL::power(x, 2);
   for(int i = 1; i < NTL::deg(poly)+1; i++)
      temp += i * poly.rep[i] * NTL::power(x, i - 1);
   
   return temp;
}
// }}}

// {{{ toP
//Converts a large integer(ZZ) to an integer modulo some prime (ZZ_p)
NTL::ZZ_p toP(NTL::ZZ & myInt)
{
   int i, nBits;
   NTL::ZZ_p temp, pow2;

   temp = 0;
   pow2 = 1;
   nBits = NTL::NumBits(myInt);
   for(i = 0; i < nBits; i++)
   {
      temp += pow2 * bit(myInt, i);
      pow2 *= 2;
   }
   if(myInt < 0)
   {
      temp *= -1;
   }
   return temp;
}
// }}}

// {{{ norm()
//Algebraic Norm(a,b) = (-b)^d * f(-a/b)
NTL::ZZ norm(int a, int b, NTL::ZZ ** poly, int d)
{
   // TODO: degree 3 specific
   // Note: Norm is specific towards f of degree 3
   NTL::ZZ temp, negOne, A, B;

   negOne = -1;
   A = a;
   B = b;
   temp = NTL::power(A, d);
   for(int i = 0; i < d; i++)
   {
      temp +=
         NTL::power(negOne, i + 3) * 
            (*(poly[i])) * NTL::power(A, i) * NTL::power(B, 3 - i);
   }
   return temp;
}
// }}}

// {{{ squareMult()
// Computes z^(limit) modulo f using the square & multiply algorithm
NTL::ZZ_pX squareMult(NTL::ZZ & limit, NTL::ZZ_pX & z, NTL::ZZ_pX & f)
{
   NTL::ZZ_pX result;

   result = 1;
   if(limit <= 0)
   {
      return result;
   }
   for(int j = NTL::NumBits(limit) - 1; j >= 0; j--)
   {
      result = (result * result) % f;
      if(bit(limit, j) == 1)
      {
         result = (z * result) % f;
      }
   }
   return result;
}
// }}}

// {{{ dfAFB()
// Returns the derivative of the polynomial as a polynomial
NTL::ZZ_pX dfAFB(NTL::ZZX &poly)
{
   // TODO: degree 3 specific
   NTL::ZZ_pX x2(2, 1), x1(1, 1);

   return ((NTL::deg(poly)+1) * x2 + toP(poly.rep[2])*2*x1 + toP(poly.rep[1]));
}
// }}}

// {{{ delta_I()
NTL::RR delta_I(NTL::ZZX &poly, int i, NTL::RR root)
{
   NTL::RR temp;

   temp = 0;
   for(int j = 0; j <= NTL::deg(poly)+1 - i - 1; j++)
      temp += NTL::MakeRR(poly.rep[i + j], 0) * NTL::power(root, j);
   
   temp += NTL::power(root, NTL::deg(poly)+1 - i);
   return temp;
}
// }}}

// {{{ factorizeRFB()
// Used to build the factorization vector for y^2 over the RFB
void factorizeRFB(int a, int b, std::vector<int> &RFB, int t, int *RFBvec, 
   NTL::ZZ & m)
{
   NTL::ZZ valZ;

   valZ = a + b * m;
   if(valZ < 0)
   {
      RFBvec[0]++;
   }

   for(int i = 0; i < t; i++)
   {
      if(valZ % RFB[i] == 0)
      {
         valZ /= RFB[i];
         RFBvec[1 + i--]++;
      }
   }
}
// }}}

// {{{ factorizeAFB()
// Used to build the factorization vector for x^2 over the AFB
void factorizeAFB(int a, int b, std::vector<int> &AFB, std::vector<int> &AFBr,
   int u, int *AFBvec, GNFS::Polynomial &polynomial)
{
   NTL::ZZ valZ;

   valZ = algebraic_norm(polynomial, a, b);
   if(valZ < 0)
   {
      AFBvec[0]++;
   }

   for(int i = 0; i < u; i++)
   {
      if((valZ % AFB[i] == 0) && ((a + b * AFBr[i]) % AFB[i] == 0))
      {
         valZ /= AFB[i];
         AFBvec[1 + i--]++;
      }
   }
}
// }}}

// {{{ normOFdf()
// Algebraic norm of the derivative of the polynomial
NTL::ZZ normOFdf(NTL::ZZX &poly)
{
   NTL::RR temp, root, r2, a2, a1, a0;
   NTL::ZZ val;

   val = poly.rep[2];
   a2 = NTL::MakeRR(val, 0);
   val = poly.rep[1];
   a1 = NTL::MakeRR(val, 0);
   val = poly.rep[0];
   a0 = NTL::MakeRR(val, 0);
   root = (-a2 + SqrRoot(abs(a2 * a2 - 3 * a1))) / 3;
   temp = NTL::power(root, 3) + a2 * NTL::power(root, 2) + a1 * root + a0;
   root = (-a2 - SqrRoot(NTL::abs(a2 * a2 - 3 * a1))) / 3;
   temp *= NTL::power(root, 3) + a2 * NTL::power(root, 2) + a1 * root + a0;
   temp *= 27;
   return NTL::RoundToZZ(temp);
}
// }}}


// {{{ square_root()
void square_root(GNFS::Polynomial &polynomial, GNFS::Target &target, Matrix &matrix,
	int pairs_needed, int nRows,  FactorBase &fb, 
	const std::vector<int> &av, const std::vector<int> &bv,
   NTL::ZZ &xZ, NTL::ZZ &yZ)
{
   using namespace std;
   using namespace NTL;

   NTL::ZZ numZ;
   NTL::ZZ valZ;
   int i;
   int j;
   bool stop = false;
   int attempt = 0;
   ZZ **piArray;                // Array of addresses of primes for Square Root
   int numPrimes, MAXP;         // Number of primes in Square Root Step
   int u = polynomial.d*target.t;


   for(int number = -1; number < nRows && !stop && attempt < 15; number++)
   {
      cout << "\tAttempt #" << ++attempt << endl;

      while(matrix.sfreeCols[++number]==0)
      {
      }

      // Store dependent (a,b) 
      for(i = 0; i < pairs_needed; i++)
      {
         matrix.sdependent[i] = false;
      }
      for(i = 0; i < nRows; i++)
      {
         if(matrix.sM[i][i]==1 && matrix.sM[i][number]==1)
         {
            matrix.sdependent[i] = true;
         }
      }
      matrix.sdependent[number] = true;


      // Calculate y
      int *RFBvec;
      RFBvec = (int *)malloc((target.t + 1) * sizeof(int));
      for(i = 0; i <= target.t; i++)
      {
         RFBvec[i] = 0;
      }
      for(i = 0; i < pairs_needed; i++)
      {
         if(matrix.sdependent[i]==1)
         {
            factorizeRFB(av[i], bv[i], fb.RFB, target.t, RFBvec, polynomial.m);
         }
      }
      yZ = 1;
      for(i = 0; i <= target.t; i++)
      {
         RFBvec[i] /= 2;
      }
      for(i = 1; i <= target.t; i++)
      {
         for(j = RFBvec[i]; j > 0; j--)
         {
            yZ = MulMod(yZ, fb.RFB[i - 1], target.n);
         }
      }
      if(RFBvec[0] % 2 == 1)
      {
         yZ *= -1;
      }


      // s(m) * f'(m)

      ZZ dF_temp = dF(polynomial.f, polynomial.m);
      std::cout << dF_temp << std::endl;

      /*
      ZZ dF_temp2 = NTL_extension::ZZX_evaluate(polynomial.f, polynomial.m);
      std::cout << dF_temp2 << std::endl;

      ZZ dF_temp3 = NTL_extension::ZZX_evaluate(NTL::diff(polynomial.f), polynomial.m);
      std::cout << dF_temp3 << std::endl;
      */

      yZ = MulMod(dF_temp, yZ, target.n);
      free(RFBvec);

      // Calculate x

      // Local Vars
      ZZ p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
      ZZ prime;

      ZZ_p::init(target.n);
      ZZ_p x;
      RR rem;

      MAXP = 10;
      piArray = (ZZ **) malloc(MAXP * sizeof(ZZ *));

      piArray[0] = &p1;
      piArray[1] = &p2;
      piArray[2] = &p3;
      piArray[3] = &p4;
      piArray[4] = &p5;
      piArray[5] = &p6;
      piArray[6] = &p7;
      piArray[7] = &p8;
      piArray[8] = &p9;
      piArray[9] = &p10;

      // Determine compatible primes

      RR product;
      RR b0, b1, b2;

      b0 = 0;
      b1 = 0;
      b2 = 0;
      for(i = 0; i < 3; i++)
      {
         product = 1;
         for(j = 0; j < pairs_needed; j++)
         {
            if(matrix.sdependent[j]==1)
            {
               product *= av[j] + NTL::to_RR(polynomial.m) * bv[j];
            }
         }
         product = SqrRoot(abs(product));
         b0 += delta_I(polynomial.f, 1, NTL::to_RR(polynomial.m))*product;
         b1 += delta_I(polynomial.f, 2, NTL::to_RR(polynomial.m))*product;
         b2 += delta_I(polynomial.f, 3, NTL::to_RR(polynomial.m))*product;
      }


      product = MakeRR(polynomial.m, 0);
      product = abs(b0) + product * abs(b1) + product * product * abs(b2);
      product *= 2.01;

      // Possible Precision Problem
      RR::SetPrecision(NumBits(TruncToZZ(product)));

      numPrimes = polynomial.d;

      numZ = NumBits(target.n);
      valZ = NumBits(TruncToZZ(product)) / numPrimes;

      while(valZ > numZ)
      {
         valZ = NumBits(TruncToZZ(product)) / ++numPrimes;
         if(numPrimes > MAXP)
         {
            valZ = 1;
            numPrimes = MAXP;
         }
      }

      numZ = numPrimes;
      b0 = 1 / MakeRR(numZ, 0);
      product = pow(product, b0);

      prime = TruncToZZ(product);

      cout << "\tBits in prime: " << NumBits(prime) << endl;

      cout << "\tFinding compatible primes" << endl;
      for(i = 0; i < numPrimes; i++)
      {
         prime = GenPrime_ZZ(NumBits(prime), 80);
         ZZ_p::init(prime);
         ZZ_pX f(3, 1), g(2, 1), h(1, 1), func;

         f = f + g * toP(polynomial.f.rep[2]) + 
				h * toP(polynomial.f.rep[1]) + toP(polynomial.f.rep[0]);

         func = (squareMult(prime, h, f) - h) % f;
         if(GCD(f, (func)) == 1)
         {
            switch (i)
            {
               case 0:
                  p1 = prime;
                  break;
               case 1:
                  p2 = prime;
                  break;
               case 2:
                  p3 = prime;
                  break;
               case 3:
                  p4 = prime;
                  break;
               case 4:
                  p5 = prime;
                  break;
               case 5:
                  p6 = prime;
                  break;
               case 6:
                  p7 = prime;
                  break;
               case 7:
                  p8 = prime;
                  break;
               case 8:
                  p9 = prime;
                  break;
               case 9:
                  p10 = prime;
                  break;
            }
         }
         else
         {
            i--;
         }
      }

      cout << "\tFinding Square Root" << endl;

      x = 0;
      rem = 0;

      for(i = 0; i < numPrimes; i++)
      {

         // Local Vars
         ZZ limit, s;
         int r;

         prime = (*(piArray[i]));
         ZZ_p::init(prime);

         ZZ_pX f(3, 1), g(2, 1), h(1, 1), z;
         ZZ_pX result, zero;

         // Setup f
         f = f + g * toP(polynomial.f.rep[2]) + h * 
				toP(polynomial.f.rep[1]) + toP(polynomial.f.rep[0]);

         limit = (prime * prime * prime - 1) / 2;
         s = limit;
         r = 1;
         while(s % 2 == 0)
         {
            s /= 2;
            r++;
         }

         result = h;
         int numTimes = 0;

         while((result + 1) != zero)
         {
            z = h + numTimes++;
            result = squareMult(limit, z, f);
         }
         z = squareMult(s, z, f);

         // Calculate delta 
         ZZ_pX delta;
         ZZ_pX gamma, deltaS;

         ZZ_p Ap, Bp;

         delta = 1;
         for(j = 0; j < pairs_needed; j++)
         {
            if(matrix.sdependent[j]==1)
            {
               Ap = av[j];
               Bp = bv[j];
               delta = (delta * (Ap + Bp * h)) % f;
            }
         }
         delta = (delta * dfAFB(polynomial.f) * 
				dfAFB(polynomial.f)) % f;
         deltaS = squareMult(s, delta, f);

         // Find gamma^2 = delta^s
         gamma = 1;
         for(j = 0; j < pow((double)2, (double)r); j++)
         {
            if((gamma * gamma) == (deltaS))
            {
               break;
            }
            gamma *= z;
         }
         limit = (s + 1) / 2;
         ZZ_pX v;

         v = squareMult(limit, delta, f) / ConstTerm(gamma);

         int *AFBvec;
         AFBvec = (int *)malloc((u + 1) * sizeof(int));
         for(j = 0; j <= u; j++)
         {
            AFBvec[j] = 0;
         }
         for(j = 0; j < pairs_needed; j++)
         {
            if(matrix.sdependent[j]==1)
            {
               factorizeAFB(av[j], bv[j],
                     fb.AFB, fb.AFBr, u, AFBvec, polynomial);
            }
         }
         for(j = 0; j <= u; j++)
         {
            AFBvec[j] /= 2;
         }
         ZZ_p normP;

         normP = 1;
         for(j = 1; j <= u; j++)
         {
            for(int k = AFBvec[j]; k > 0; k--)
            {
               normP *= fb.AFB[j - 1];
            }
         }
         valZ = normOFdf(polynomial.f);
         normP *= toP(valZ);
         limit = (prime * prime * prime - 1) / (prime - 1);
         if(ConstTerm(squareMult(limit, v, f)) != normP)
         {
            v *= -1;
         }
         free(AFBvec);

         ZZ_p xi, ai, Pi;

         xi = toP(polynomial.m) * toP(polynomial.m) * coeff(v, 2) 
				+ toP(polynomial.m) * coeff(v, 1) + coeff(v, 0);
         Pi = 1;
         for(j = 0; j < numPrimes; j++)
         {
            if(j != i)
            {
               Pi *= toP((*(piArray[j])));
            }
         }
         ai = 1 / Pi;
         ZZ_p::init(target.n);
         Pi = 1;
         for(j = 0; j < numPrimes; j++)
         {
            if(j != i)
            {
               Pi *= toP((*(piArray[j])));
            }
         }

         x += ai * xi * Pi;

         RR::SetPrecision(2 * NumBits(prime));

         RR aiR, xiR, piR;

         // Possible Precision Problem
         aiR = MakeRR(rep(ai), 0);
         xiR = MakeRR(rep(xi), 0);
         piR = MakeRR((*(piArray[i])), 0);

         rem += aiR * xiR / piR;

      }

      ZZ_p::init(target.n);
      ZZ_p P;

      P = 1;
      for(j = 0; j < numPrimes; j++)
      {
         P *= toP((*(piArray[j])));
      }
      valZ = RoundToZZ(rem);
      x -= toP(valZ) * P;

      stop = false;
      if(GCD(rep(x) - yZ, target.n) != 1 && GCD(rep(x)-yZ,target.n)!=target.n)
      {
         cout << "\n\ty = " << yZ;
         cout << "\n\tx = " << x << endl;
         xZ = rep(x);
         stop = true;
      }
      else if(GCD(rep(x)+yZ,target.n) != 1 && GCD(rep(x)+yZ,target.n)!=target.n)
      {
         cout << "\n\ty = " << yZ;
         cout << "\n\tx = " << x << endl;
         xZ = rep(x);
         stop = true;
      }

      free(piArray);
      if(!stop)
      {
         cout << "\tNo Factor Found" << endl;
         //exit(0);
      }

   }
   // Ends loop over different (a,b) pairs from
}
// }}}



