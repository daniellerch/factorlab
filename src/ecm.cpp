/*
   Copyright 2015 Daniel Lerch Hostalot.

   This file is part of FactorLib.

   FactorLib is free software; you can redistribute it and/or modify                
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



#include <factorlib.hpp>

#include <ec/EC_p.h>

using namespace NTL;

// {{{ ECM_parameters()
// get optimal parameters                                                        
void ECM_parameters(long& B1, long& B2, double& prob, long& D,                   
          long pbits, long nbits) {                                              
  // smoothness bounds                                                           
  double logp = pbits*M_LN2;                                                     
  double b1 = 1.25*exp(0.77*sqrt(logp)*sqrt(log(logp)));                         
  B1 = b1<NTL_SP_BOUND ? (long)b1 : NTL_SP_BOUND;                                
  double b2 = 70*exp(0.78*sqrt(logp)*sqrt(log(logp)));                           
  B2 = b2<NTL_SP_BOUND ? (long)b2 : NTL_SP_BOUND;                                
                                                                                 
  // memory use in stage 2 is 3*D*nbits/8                                        
  D = (long)sqrt((B2-B1)/2.0);                                                   
  if (2*D>=B1) D = (B1-1)/2;                                                     
  if (D<2) D=2;                                                                  
                                                                                 
  // all curve orders are divisible by 12                                        
  logp-=log(12.0);                                                               
                                                                                 
  // probability that curve order is B1-smooth                                   
  double logB1 = log((double)B1);                                                
  double u = logp/logB1;                                                         
  prob = std::pow(u,-u);                                                         
  // add in probability that order has one larger factor (up to B2)              
  // numerical integration: time is O(log(B2/B1))                                
  long min_x = B1;                                                               
  double min_y = std::pow(u-1,1-u)/logB1/B1;                                     
  while (min_x<B2) {                                                             
    long max_x = 2*min_x;                                                        
    if (max_x>B2 || min_x>max_x) max_x=B2;                                       
    double logx = log((double)max_x);                                            
    double v = logx/logB1;                                                       
    double max_y = std::pow(u-v,v-u)/logx/max_x;                                 
    prob += (max_x-min_x)*(min_y+max_y)/2;                                       
    min_x = max_x;                                                               
    min_y = max_y;                                                               
  }                                                                              
  if (prob>1) prob=1;                                                            
} 
// }}}

// returns a non-trivial factor q of ZZ_p::modulus(), or 1 otherwise             
void ECM_stage_one(ZZ& q, EC_p& Q, PrimeSeq& seq, long bound) {
  long sbound = (long)sqrt((double)bound);                                       
  seq.reset(0);                                                                  
  long p = seq.next();                                                           
  for (; p<=sbound; p=seq.next()) {                                              
    long pp,t=p;                                                                 
    do { pp=t; t*=p; } while (t>pp && t<=bound); // we might overflow t here     
    mul(Q,Q,pp);                                                                 
  }                                                                              
  for (; p<=bound; p=seq.next())                                                 
    mul(Q,Q,p);                                                                  
  if (!IsZero(Q))                                                                
    GCD(q,ZZ_p::modulus(),rep(Q.Z));                                             
  else                                                                           
    set(q);                                                                      
}  

// returns a non-trivial factor q of ZZ_p::modulus(), or 1 otherwise             
void ECM_stage_two(ZZ& q, EC_p& Q, PrimeSeq& seq, long bound, long D) {          
  long B1 = seq.next();                                                          
  if (B1<=2*D) {                                                                 
    // no primes to work with                                                    
    set(q);                                                                      
    return;                                                                      
  }                                                                              
                                                                                 
  EC_p R;                                                                        
  mul(R,B1,Q);                                                                   
  // check R for divisor                                                         
  if (IsZero(R)) {                                                               
    set(q);                                                                      
    return;                                                                      
  }                                                                              
  GCD(q,ZZ_p::modulus(),rep(R.Z));                                               
  if (!IsOne(q))                                                                 
    return;                                                                      
  EC_p T;                                                                        
  mul(T,B1-2*D,Q);                                                               
                                                                                 
  // compute point multiples S[d]=2*d*Q                                          
  EC_p S[D+1];                                                                   
  S[0]=Q;                                                                        
  doub(S[1],Q);                                                                  
  doub(S[2],S[1]);                                                               
  for (long d=3; d<=D; ++d)                                                      
    addh(S[d],S[d-1],S[1],S[d-2]);                                               
  ZZ_p beta[D+1];                                                                
  for (long d=0; d<=D; ++d)                                                      
    mul(beta[d],S[d].X,S[d].Z);                                                  
       
  ZZ_p g,t,t2;                                                                   
  set(g);                                                                        
  ZZ_p alpha;                                                                    
  long r=B1;                                                                     
  long p=seq.next();                                                             
  do {                                                                           
    mul(alpha,R.X,R.Z);                                                          
    do {                                                                         
      long delta = (p-r)/2;                                                      
      if (delta>D) break;                                                        
      //g *= (R.X-S[delta].X) * (R.Z+S[delta].Z) - alpha + beta[delta];          
      sub(t,R.X,S[delta].X);                                                     
      add(t2,R.Z,S[delta].Z);                                                    
      t*=t2;                                                                     
      t-=alpha;                                                                  
      t+=beta[delta];                                                            
      g*=t;                                                                      
      // next prime                                                              
      p = seq.next();                                                            
      if (p==0) {                                                                
   // ran out of primes (should never happen)                                    
   p=NTL_MAX_LONG;                                                               
   break;                                                                        
      }                                                                          
    } while (p<=bound);                                                          
    if (p>bound)                                                                 
      break;                                                                     
    addh(T,R,S[D],T);                                                            
    swap(R,T);                                                                   
    r+=2*D;                                                                      
  } while (true);                                                                
  if (!IsZero(g))                                                                
    GCD(q,ZZ_p::modulus(),rep(g));                                               
  else                                                                           
    set(q);                                                                      
}

void ECM_random_curve(EC_pCurve& curve, ZZ_p& X, ZZ_p& Z) {                      
  ZZ sigma;                                                                      
  RandomBnd(sigma,ZZ_p::modulus()-6);                                            
  sigma+=6;                                                                      
  ZZ_p u,v;                                                                      
  u = sqr(to_ZZ_p(sigma))-5;                                                     
  v = 4*to_ZZ_p(sigma);                                                          
  ZZ_p C,Cd;                                                                     
  C = (v-u)*sqr(v-u)*(3*u+v);                                                    
  Cd = 4*u*sqr(u)*v;                                                             
  // make sure Cd is invertible                                                  
  ZZ Cinv;                                                                       
  if (InvModStatus(Cinv,rep(Cd),ZZ_p::modulus())!=0) {                           
    conv(X,Cinv);                                                                
    clear(Z);                                                                    
    return;                                                                      
  }                                                                              
  C*=to_ZZ_p(Cinv);                                                              
  C-=2;                                                                          
                                                                                 
  // random curve                                                                
  ZZ_pX f;                                                                       
  SetCoeff(f,3);                                                                 
  SetCoeff(f,2,C);                                                               
  SetCoeff(f,1);                                                                 
  conv(curve,f);                                                                 
  curve.SetRepresentation(curve.MONTGOMERY);                                     
                                                                                 
  // initial point                                                               
  mul(X,u,sqr(u));                                                               
  mul(Z,v,sqr(v));                                                               
}  
void ECM_one_curve(ZZ& q, PrimeSeq& seq, long B1, long B2, long D) {             
  EC_pBak bak;                                                                   
  bak.save();                                                                    
                                                                                 
  // random curve and point                                                      
  EC_pCurve curve;                                                               
  ZZ_p X,Z;                                                                      
  ECM_random_curve(curve,X,Z);                                                   
  if (IsZero(Z)) {                                                               
    q=rep(X);                                                                    
    return;                                                                      
  }                                                                              
  EC_p::init(curve);                                                             
                                                                                 
  // initial point                                                               
  EC_p Q;                                                                        
  Q.X = X;                                                                       
  Q.Z = Z;                                                                       
                                                                                 
  // attempt to factor                                                           
  ECM_stage_one(q,Q,seq,B1);                                                     
  if (IsOne(q)&&!IsZero(Q))                                                      
    ECM_stage_two(q,Q,seq,B2,D);                                                 
}  




// {{{ FL::ecm()
NTL::ZZ FL::ecm(NTL::ZZ &n)
{
   NTL::ZZ B = NTL::to_ZZ(10000000);

   long B1,B2,D;                                                                  
   double prob;                                                                   
   ECM_parameters(B1, B2, prob, D, NTL::NumBits(B), NTL::NumBits(n)); 

   NTL::ZZ_p::init(n);

   NTL::PrimeSeq seq;                                                                  

   for (long i=0; i<100; ++i) 
   {
      ECM_one_curve(n, seq, B1, B2, D);                                                
      if (!NTL::IsOne(n)) 
      {  
         std::cout << "factor:" << n << std::endl;
         return n;
      }                                                                            
   } 
                                                                             
   // TODO

   std::cout << "Warning: not implemented!" << std::endl;
   return n;   
}
// }}}


