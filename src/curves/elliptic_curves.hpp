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




#ifndef __ELLIPTIC_CURVES_HPP__
#define __ELLIPTIC_CURVES_HPP__

#include <ostream>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>



namespace EC
{
 
   /* Montgomery Curves and Points */

	// {{{ Curve
   // g*y^2 = x^3 + Cx^2 + Ax + B
	class Curve 
	{ 
		public:
			NTL::ZZ_pX f;

		public:
			Curve() {}
			~Curve() {}

         void set(const NTL::ZZ_p &A, const NTL::ZZ_p &B, const NTL::ZZ_p &C)
         {
            NTL::ZZ_pX newf;
            NTL::SetCoeff(newf, 3); // monic
            NTL::SetCoeff(newf, 2, C);
            NTL::SetCoeff(newf, 1, A);
            NTL::SetCoeff(newf, 0, B);
            f=newf;
         }

         void get(NTL::ZZ_p &A, NTL::ZZ_p &B, NTL::ZZ_p &C) const
         {
            NTL::GetCoeff(C, f, 2);
            NTL::GetCoeff(A, f, 1);
            NTL::GetCoeff(B, f, 0);
         }


			const NTL::ZZ_pX& get_poly() const 
         {
				return f;
			}

			Curve& operator=(const Curve& c) 
         {
				f=c.f;
				return *this;
			}
	};
	// }}}

	// {{{ Point
	class Point 
	{
      public:
         Curve curve;
			NTL::ZZ_p X, Z;

		public:

			Point() {}
			Point(const Curve &c) { curve=c; }
			~Point() {}

         Curve get_curve()
         {
            return curve;
         }

         Point& operator=(const Point& a)
         {
            NTL::ZZ_p A, B, C;
            a.curve.get(A, B, C);
            curve.set(A, B, C);
            X=a.X; 
            Z=a.Z;
            return *this;
         }
	};
	// }}}


   Point operator*(const Point& a, const NTL::ZZ& b);
   std::ostream& operator<<(std::ostream& s, const Curve& curve);
   bool operator!=(const Point& a, const Point& b);
   bool operator==(const Point& a, const Point& b);
   std::ostream& operator<<(std::ostream& s, const Point& x);


	/* Arithmetic */

   // {{{ doubleh()
   // g*y2 = x3 + Cx2 + Ax + B
   // 7.7 : Prime Numbers. Crandall & Pomerance.
   inline void doubleh(Point& x, const Point& a) 
   {
      if(a.Z==0) 
      {
         // Infinity
         x.X=0;
         x.Z=0;
         return;
      }

      NTL::ZZ_p X, Z;

      NTL::ZZ_p A, B, C;
      a.curve.get(A, B, C);
      x.curve.set(A, B, C);

      NTL::ZZ_p X2, Z3;
      NTL::sqr(X2, a.X);
      Z3 = a.Z*a.Z*a.Z;

      // X+ = (X1^2-AZ1^2)^2 -4B(2X1+CZ1)Z1^3
      X = sqr(X2-A*sqr(a.Z)) - 4*B*(2*a.X+C*a.Z)*Z3;

      // Z+ = 4Z1(X1^3+CX1^2 Z1+AX1Z1^2+BZ1^3)
      Z = 4*a.Z*(a.X*X2 + C*X2*a.Z + A*a.X*sqr(a.Z) + B*Z3);

      x.X = X;
      x.Z = Z;
   }

   inline Point doubleh(const Point& a) 
   {
      Point x(a.curve);  
      doubleh(x,a);  
      return x;
   }

   // }}}

     // {{{ addh()
   // g*y2 = x3 + Cx2 + Ax + B
   // 7.6 : Prime Numbers. Crandall & Pomerance.
   inline void addh(Point& x, const Point& a, const Point& b, const Point& c) 
   {
      NTL::ZZ_p X, Z, A, B, C;

      a.curve.get(A, B, C);
      x.curve.set(A, B, C);

      NTL::ZZ_p Z1Z2 = a.Z*b.Z;

      // X+ = Z_ ((X1X2-AZ1Z2)^2 - 4B(X1Z2+X2Z1+CZ1Z2)Z1Z2)
      X = c.Z * (NTL::sqr(a.X*b.X - A*Z1Z2) - 4*B*(a.X*b.Z + b.X*a.Z + C*Z1Z2)*Z1Z2);

      // Z+ = X_ (X1Z2-X2Z1)^2         
      Z = c.X * NTL::sqr(a.X*b.Z - b.X*a.Z);

      x.X=X;
      x.Z=Z;
   }
   // }}}

   // {{{ mul()
   // Algorithm 7.2.7 : Prime Numbers. Crandall & Pomerance.
   inline void mul(Point& x, const Point& P, const NTL::ZZ& k) 
   {
      if(k==0)
      {
         // Infinity
         x.X=0;
         x.Z=0;
         return;
      }

      if(k==1)
      {
         // Return the original point
         x=P;
         return;
      }

      if(k==2)
      {
         doubleh(x, P);
         return;
      }

      // Begin  adding/doubling ladder
      Point t(P.curve);
      Point u = P;
      doubleh(t, P);

      // Loop over bits of n, starting with next-to-highest
      for(long j=NumBits(k)-2; j>0; --j) 
      {
         if(NTL::bit(k, j) == 1) 
         {
            addh(u, t, u, P);
            doubleh(t, t);
         }
         else    
         {
            addh(t, u, t, P);
            doubleh(u,u);
			}
		}

      // Final calculation
		if(NTL::IsOdd(k))
			addh(x, u, t, P);
		else
			doubleh(x, u);
	}

	inline void mul(Point& x, const Point& P, long k) 
	{
		mul(x, P, NTL::to_ZZ(k));
	}

   inline void mul(Point& x, const NTL::ZZ& a, const Point& b) 
   {
      mul(x, b, a);
   }

   inline void mul(Point& x, long a, const Point& b) 
   {
      mul(x, b, a);
   }
	// }}}


	/* Operators */

   // {{{ operator * 
   inline Point operator*(const Point& a, const NTL::ZZ& b) 
   {
      Point x(a.curve);  
      mul(x, a, b);  
      return x;
   }

   inline Point operator*(const Point& a, long b) 
   {
      Point x(a.curve);  
      mul(x, a, b);
      return x;
   }

   inline Point operator*(const NTL::ZZ& a, const Point& b) 
   {
      Point x(b.curve);
      mul(x, a, b);
      return x;
   }

   inline Point operator*(long a, const Point& b) 
   {
      Point x(b.curve);
      mul(x, a, b);
      return x;
   }
   // }}}

   // {{{ operator << 
   inline std::ostream& operator<<(std::ostream& s, const Curve& curve) 
   {
      return s << curve.get_poly();
   }
   // }}}

   // {{{ operator != 
   inline bool operator!=(const Point& a, const Point& b) 
   {
      return !(a==b);
   }
   // }}}

   // {{{ operator == 
   inline bool operator==(const Point& a, const Point& b) 
   {
      if(IsZero(a.Z))
         return IsZero(b.Z);

      if(IsZero(b.Z))
         return false;

      if(a.Z==b.Z) 
         return a.X==b.X;

      return a.X*b.Z==b.X*a.Z;
   }
   // }}}

   // {{{ operator << 
   inline std::ostream& operator<<(std::ostream& s, const Point& x) 
   {
      return s << "[" << x.X << " : " << x.Z <<"]";
   }
	// }}}

}

#endif
