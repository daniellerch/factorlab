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

// {{{ get_k()
static NTL::ZZ get_k(long B1)
{
   NTL::PrimeSeq seq;
   seq.reset(0);
   NTL::ZZ r, pa;
   r=1;

   for(long p=seq.next(); p<=B1; p=seq.next())
   {
      // TODO: This can be cached for the next curve
      pa = p;
      while(pa<=B1)
         pa *= p;
      pa/=p;

      r*=pa;
   }

   return r;
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

// {{{ point_order()
// Brute force point order
long point_order(const EC::Point &P)
{
   NTL::ZZ x;
   long order=1;
   EC::Point R=P;
   for(; ; ++order)
   {
      R=order*P;
      //std::cout << order << " | " << "[" << R.X << ":" << R.Z << "]" << std::endl;

      if(R.Z!=0)
      {
         NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(R.Z));
         if(1<x && x< NTL::ZZ_p::modulus())
            return order;         
      }
   }
 
   return 0;
}
// }}}

// {{{ ecm_rho()
NTL::ZZ ecm_rho(const EC::Point &P, const NTL::ZZ &k, long max_it)
{
   NTL::ZZ m1, m2, x;
   EC::Point U=P, Uc=P;
   EC::Point V=P, Vc=P;
   EC::Point Q;
   m1=1;
   m2=1;
   x=1;

   //for(long l=0; /*l<max_it*/; l++)
   for(long l=0; l<max_it; l++)
   {
      Uc=U;
      Vc=V;
      U=k*U;
      V=k*k*V;

      NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(U.Z));
      if(1<x && x<NTL::ZZ_p::modulus())
      {
         std::cout << "U casual solution: " << l << std::endl;
         return x;
      }

      NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(V.Z));
      if(1<x && x<NTL::ZZ_p::modulus())
      {
         std::cout << "V casual solution: " << l << std::endl;
         return x;
      }




      //m1*=k;
      //m2*=k*k;

      //m1=NTL::MulMod(k, m1, NTL::ZZ_p::modulus());
      //m2=NTL::MulMod(k*k, m2, NTL::ZZ_p::modulus());

      /*
      m1=NTL::power(k, l);
      m2=NTL::power(k, 2*l);
      Q = (m1-m2)*P;
      std::cout << "mul-1:" << Q << std::endl;

      Q = (-m1+m2)*P;
      std::cout << "mul-2:" << Q << std::endl;

      Q = (m1+m2)*P;
      std::cout << "mul+ :" << Q << std::endl;

      EC::Point T, R;
      T = m1*P;
      R = m2*P;
      addh(Q, T, R, P);
      std::cout << "addh1:" << Q << std::endl;

      addh(Q, Uc, Vc, U);
      std::cout << "addh2:" << Q << std::endl;
      */

      addh(Q, Uc, Vc, P);

      NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(Q.Z));
      if(1<x && x<NTL::ZZ_p::modulus())
      {
         std::cout << "cycle (1): " << l << std::endl;
         return x;
      }

      continue;

      if(U==V)
      {
         std::cout << "cycle: " << l << std::endl;
         Q = (m1-m2)*P;
         NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(Q.Z));
         if(1<x && x<NTL::ZZ_p::modulus())
            return x;

         Q = (m1+m2)*P;
         NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(Q.Z));
         if(1<x && x<NTL::ZZ_p::modulus())
            return x;

         Q = (-m1+m2)*P;
         NTL::GCD(x, NTL::ZZ_p::modulus(), NTL::rep(Q.Z));
         if(1<x && x<NTL::ZZ_p::modulus())
            return x;
      }
   }

   return x;
}
// }}}

// {{{ ecrho()
NTL::ZZ ecrho(NTL::ZZ &n, long B1, long B2)
{
   NTL::ZZ B = NTL::SqrRoot(n);
   long curves=0, D;
   NTL::ZZ_p::init(n);
   NTL::SetSeed(NTL::to_ZZ(0));

   // TODO: Improve parameter selection
   if(B1==0) B1=2000;

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
      //std::cout << "Q order: " <<  point_order(Q) << std::endl;

      //NTL::ZZ k = get_k(B1)*4987;
      NTL::ZZ k = get_k(B1);
      std::cout << "Using k: " << k << std::endl;
      NTL::ZZ factor = ecm_rho(Q, k, B1*10);
         
      if(1<factor && factor<NTL::ZZ_p::modulus())
         return factor;
   }

   return n;
}
// }}}



