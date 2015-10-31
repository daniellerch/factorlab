/*
   Copyright 2015 Daniel Lerch Hostalot.

   This file is part of factorlab.

   factorlab is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   FLINT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FLINT; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */


#include <curves/elliptic_curves.hpp>
#include <ecm/ecm.hpp>
#include <NTL/RR.h>

// {{{ stage_one()
// Algorithm 7.4.4 (3). Prime Numbers. Crandall & Pomerance.
static void stage_one(NTL::ZZ& q, EC::Point& Q, long B1)
{
   NTL::PrimeSeq seq;
   seq.reset(0);

   for(long p=seq.next(); p<=B1; p=seq.next())
   {
      // TODO: This can be cached for the next curve
      long pa = p;
      while(pa<=B1)
         pa *= p;
      pa/=p;

      Q = pa * Q;
   }

   if (Q.Z!=0)
      NTL::GCD(q, NTL::ZZ_p::modulus(), NTL::rep(Q.Z));
   else
      q=1;
}
// }}}

// {{{ stage_two()
// Algorithm 7.4.4 (4). Prime Numbers. Crandall & Pomerance.
static void stage_two(NTL::ZZ& x, EC::Point& Q, long B1, long B2, long D)
{
   EC::Point S[D+1], T, R, Tmp;
   NTL::ZZ_p alpha, beta[D+1], g;

   doubleh(S[1], Q);
   doubleh(S[2], S[1]);

   // Compute S[d] = 2d*Q
   for (long d=0; d<=D; ++d)
   {
      if(d>2)
         addh(S[d], S[d-1], S[1], S[d-2]);

      // Store the XZ products
      beta[d] = S[d].X * S[d].Z;
   }

   long delta;
   long B=B1-1;
   g = 1;
   T = (B-2*D)*Q;
   R = B*Q;

   for(long r=B; r<B2; r=r+2*D)
   {
      alpha = R.X * R.Z;
      // TODO: Slow prime generation. Use cache.
      for(long q=r+2; q<=r+2*D; q=NTL::NextPrime(q+1))
      {
         delta = (q-r)/2; // Distance to next prime
         g *= (R.X-S[delta].X) * (R.Z+S[delta].Z) - alpha + beta[delta];
      }
      
      Tmp=R;
      addh(R, R, S[D], T);
      T=Tmp;
   }
      
   NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(g));
}
// }}}

// {{{ choose_random_curve_and_get_point()
// Algorithm 7.4.4 (2) : Theorem 7.4.3. Prime Numbers. Crandall & Pomerance.
// Given a curve governed by y^2 = x^3 + C*sigma*x^2 + x, there exists a point
// with x-coordinate u^3 v^(-3) and the order is divisible by 12.
static EC::Point choose_random_curve_and_get_point(EC::Curve& curve)
{
   NTL::ZZ_p X, Z;
   NTL::ZZ_p u,v;
   NTL::ZZ_p C, Cd;
   NTL::ZZ Cinv;

   for(;;)
   {
      // Random sigma in [6, n-1]
      NTL::ZZ_p sigma;
      NTL::random(sigma);
      sigma+=6;

      u = NTL::sqr(sigma)-5;
      v = 4*sigma;
      C = NTL::power(v-u, 3)*(3*u+v);
      Cd = 4*u*NTL::sqr(u)*v;

      if(NTL::InvModStatus(Cinv, NTL::rep(Cd), NTL::ZZ_p::modulus())!=0) 
      {
         NTL::conv(X, Cinv);
         X=0;
         continue;
      }
      break;
   }

   C = C * NTL::to_ZZ_p(Cinv) - 2;
      
   // g*y^2 = x^3 + Cx^2 + Ax + B
   // set(A, B, C)
   curve.set(NTL::to_ZZ_p(1), NTL::to_ZZ_p(0), C);


   //std::cout << u << "|" << v << "|" << C << std::endl;

   // Initial point
   EC::Point Q(curve);
   Q.X = NTL::power(u, 3);
   Q.Z = NTL::power(v, 3);

   return Q;
}
// }}}

// {{{ ecm()
NTL::ZZ ecm(NTL::ZZ &n, long B1, long B2)
{
   NTL::ZZ B = NTL::SqrRoot(n);
   long curves=0, D;
   NTL::ZZ_p::init(n);
   NTL::SetSeed(NTL::to_ZZ(0));

   // TODO: Improve parameter selection
   if(B1==0) B1=2000;
   if(B2==0) B2=100*B1;

   D=(long)sqrt((B2-B1)/2.0);                                                   
   if (2*D>=B1) D = (B1-1)/2;                                                     
   if (D<2) D=2;  


   for (;;curves++)
   {
		NTL::ZZ x = n;

      EC::Curve curve;
      EC::Point Q = choose_random_curve_and_get_point(curve);
      if (IsZero(Q.Z))
         continue;
      std::cout << curves << ": Curve: " << curve << std::endl;

      NTL::ZZ m;
      for(int j=0; j<1; j++)
      {  // random points per curve
         std::cout << j << ": Point Q: " << Q << std::endl;   

         stage_one(x, Q, B1);
         if(1<x && x<n) 
         {  
            std::cout << "Factor stage one" << std::endl;
            return x;
         }

         stage_two(x, Q, B1, B2, D);
         if(1<x && x<n)
         {
            std::cout << "Factor stage two" << std::endl;
            return x;
         }

         // Chose another random point in the curve
         NTL::RandomBnd(m, NTL::ZZ_p::modulus());
         Q = m*Q;

         if(Q.Z==0)
            break;
      }
   }

   return n;
}
// }}}



