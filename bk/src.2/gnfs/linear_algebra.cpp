
#include <sieve.hpp>
#include <matrix/block_lanczos.hpp>
#include <matrix/smatrix.h>
#include <matrix/smat_long.h>

#include <linear_algebra.hpp>

#include <NTL/GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>


// {{{ get_freecols()
NTL::svec_long get_freecols(NTL::smat_long& M) 
{
   NTL::svec_long freecols;
   freecols.SetLength(M.NumCols());

   for(int i=0, j=0; (i<M.NumRows())&&(j<M.NumCols()); i++, j++) 
   {
      if(M[i][j]==0) 
      {
         int k=i+1;
         while((k<M.NumRows())&&(M[k][j]==0)) k++;

         if(k>=M.NumRows()) 
         {
            freecols[j]=1;
            i--;
         }
      }
   }  

   return freecols;
}  
// }}}

// {{{ gaussian_elimination()
void gaussian_elimination(NTL::smat_long &M)
{
   long i, j, k;
   bool stop;

   // Solve Mx = 0
   //NTL::GF2 temp;
   long temp;

   for(i=0; i < M.NumCols()-1; i++)
   {
      stop = false;
      j = i;
      if(M[i][j]==0)
      {
         // find k so that M[k,j]==1
         k = i;
         while(++k < M.NumCols()-1 && M[k][j]==0)
         {
         }
         // swap row i with row k
         if(k >= M.NumCols()-1)
         {
            stop = true;
         }
         for(j = 0; j <= M.NumRows()-1 && !stop; j++)
         {
            temp = M[i][j];
            M[i][j] = M[k][j];
            M[k][j] = temp;
         }

      }
      for(k = 0; k < M.NumCols()-1 && !stop; k++)
      {
         if(M[k][i]==1 && k != i)
         {
            for(j = 0; j <= M.NumRows()-1; j++)
            {
               //M[k][j] =  
               //((NTL::IsOne(M[k][j]) ||NTL::IsOne(M[i][j])) &&
               //((NTL::IsZero(M[k][j])||NTL::IsZero(M[i][j]))));
               M[k][j] = (M[k][j]!=0 || M[i][j]!=0) &&
                         (M[k][j]==0 || M[i][j]==0);
            }
         }
      }
   }

}
// }}}


// {{{ Legendre()
//1 if y is a quadratic residue modulo z, -1 otherwise
int Legendre(NTL::ZZ & y, NTL::ZZ & z)
{
   NTL::ZZ x, yP;

   yP = y % z;
   for(x = 0; x < z; x++)
   {
      if(((x * x) - yP) % z == 0)
      {
         return 1;
      }
   }
   return -1;
}
// }}}

// {{{ linear_algebra()
void linear_algebra(
   GNFS::Polynomial &polynomial, 
   GNFS::Target &target, 
	FactorBase &fb, 
   Matrix &matrix, 
   const std::vector<int> &av, 
   const std::vector<int> &bv)
{

   NTL::ZZ aZ;
   NTL::ZZ bZ;
   NTL::ZZ pZ;
   NTL::ZZ valZ;
   NTL::ZZ numZ;
   int i;
   int j;
   int k;
   int u = polynomial.d*target.t;


   // Initialize sM
   for(j = 0; j <= matrix.sM.NumCols()-1; j++)
   {
      // Initialize row
      for(k = 0; k < matrix.sM.NumCols()-1; k++)
         matrix.sM[k][j] = 0;
      
      // Set the first column
      aZ = av[j];
      bZ = bv[j];
      valZ = aZ + bZ * polynomial.m;
      if(valZ < 0)
      {
         valZ *= -1;
         matrix.sM[0][j] = 1;
      }

      // Set a RFB row 
      i = 0;
      while(i < target.t && valZ != 1)
      {
         pZ = fb.RFB[i];
         if(valZ % pZ == 0)
         {
            if(matrix.sM[1+i][j]==0) matrix.sM[1+i][j]=1;
            else matrix.sM[1+i][j]=0;
            valZ = valZ / pZ;
         } 
         else i++;
      }

      // Set a AFB row 
      valZ = algebraic_norm(polynomial, av[j], bv[j]);
      if(valZ < 0)
         valZ *= -1;
      
      i = 0;
      while(i<u && valZ!=1)
      {
         pZ = fb.AFB[i];
         if(valZ % pZ == 0)
         {
            numZ = fb.AFBr[i];
            //while((aZ + bZ * numZ) % pZ != 0 && i<target.t) // TODO
            while((aZ + bZ * numZ) % pZ != 0) // TODO
            {
               pZ = fb.AFB[++i];
               numZ = fb.AFBr[i];
            }

            if(matrix.sM[1+target.t+i][j]==0) matrix.sM[1+target.t+i][j]=1;
            else matrix.sM[1+target.t+i][j]=0;

            valZ = valZ / pZ;
         }
         else i++;
         
      }

      // Set a QCB row
      for(i=0; i<target.digits; i++)
      {
         numZ = fb.QCB[i];
         valZ = fb.QCBs[i];
         valZ = aZ + bZ * valZ;
         if(Legendre(valZ, numZ) != 1)
         {
            matrix.sM[1+target.t+u+i][j]=1;
         
         }
      }
   }
   
   std::cout << "\tSize: " << matrix.sM.NumRows() << "x" << matrix.sM.NumCols() << std::endl;

   gaussian_elimination(matrix.sM);
   matrix.sfreeCols = get_freecols(matrix.sM);

}
// }}}



