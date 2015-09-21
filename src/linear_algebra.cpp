
#include <sieve.hpp>
#include <matrix/block_lanczos.hpp>
#include <matrix/smatrix.h>
#include <matrix/smat_long.h>

#include <linear_algebra.hpp>

#include <NTL/GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>

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
void linear_algebra(GNFS::Polynomial &polynomial, GNFS::Target &target, 
	FactorBase &fb, Matrix &matrix, 
    const std::vector<int> &av, const std::vector<int> &bv)
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

//#define LANCZOS 
#ifdef LANCZOS
   /*
   uint64_t *nullrows;
   uint32_t extra_rels = 0; // number of opportunities to factor n 
   uint32_t ncols, nrows;

   ncols = matrix.col;
   nrows = matrix.row;

   printf("=>%u\n", nrows);
   la_col_t *M = create_matrix(nrows);
   */
   
   NTL::smat_long Aorig;
   Aorig.SetDims(matrix.row, matrix.col);
#endif

   // Initialize sM
   for(j = 0; j <= matrix.row-1; j++)
   {
      // Initialize row
      for(k = 0; k < matrix.col-1; k++)
         matrix.sM[k][j] = 0;
      
      // Set the first column
      aZ = av[j];
      bZ = bv[j];
      valZ = aZ + bZ * polynomial.m;
      if(valZ < 0)
      {
         valZ *= -1;
         matrix.sM[0][j] = 1;
      
#ifdef LANCZOS      
         //insert_col_entry(&M[0], j);
         Aorig[0][j]=1;
#endif
      }

      // Set a RFB row 
      i = 0;
      while(i < target.t && valZ != 1)
      {
         pZ = fb.RFB[i];
         if(valZ % pZ == 0)
         {
#ifdef LANCZOS      
            //if(matrix.sM[1+i][j]==0) insert_col_entry(&M[1+i], j);
            if(matrix.sM[1+i][j]==0) Aorig[1+i][j]=1;
#endif
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

#ifdef LANCZOS      
            //if(matrix.sM[1+target.t+i][j]==0) insert_col_entry(&M[1+target.t+i], j);
            if(matrix.sM[1+target.t+i][j]==0) Aorig[1+target.t+i][j]=1;
#endif
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
         
#ifdef LANCZOS      
            //insert_col_entry(&M[1+target.t+u+i], j);
            Aorig[1+target.t+u+i][j]=1;
#endif
         }
      }
   }
   
   std::cout << "\tSize: " << matrix.sM.NumRows() << "x" << matrix.sM.NumCols() << std::endl;

#ifdef LANCZOS      
   std::cout << "OLD MATRIX" << std::endl;
   std::cout << matrix.sM << std::endl;

   std::cout << "LANCZOS C++ MATRIX" << std::endl;
   for(unsigned int i=0; i<matrix.row; i++)
   {
      for(unsigned int j=0; j<matrix.col; j++)
         std::cout << Aorig[i][j] << " ";
      std::cout << std::endl;
   }

   /*
   for(unsigned int i=0; i<ncols; i++)
   {
      for(unsigned int j=0; j<nrows; j++)
      {
         int found=0;
         for(unsigned int w=0; w<M[i].weight; w++)
            if(M[i].data[w]==j) { found=1; break; }
           
         if(found) std::cout << "1 ";
         else std::cout << "0 ";
      }
      std::cout << std::endl;
   }
   */

#endif

   /*
   uint64_t *nullrows;
   uint32_t extra_rels = 64; // number of opportunities to factor n 
   uint32_t ncols, nrows;

   ncols = matrix.col + extra_rels;
   nrows = matrix.row;

   ncols = 10 + extra_rels;
   nrows = 100;

   la_col_t *M = create_matrix(nrows);

   for(uint32_t i=0; i<nrows; i++)
   {
      insert_col_entry(&M[i], rand()%ncols);
      insert_col_entry(&M[i], rand()%ncols);
      insert_col_entry(&M[i], rand()%ncols);
   }
   */


#ifdef LANCZOS      

   NTL::ZZ p = NTL::to_ZZ(2);
   NTL::ZZ_p::init(p);
   NTL::vec_long cols;
   NTL::smat_long A;                                                                   
   NTL::vec_ZZ_p y; 
   NTL::vec_ZZ_p yorig;                                                                
   yorig.SetLength(Aorig.NumRows()); 
   //for (long i=0; i<Aorig.NumRows(); ++i)
   //   yorig[i]=0;

   // Find smaller Aorig*x=yorig
   SGauss(A,y,cols,Aorig,yorig); 

   std::cout << "LANCZOS C++ MATRIX (s-gauss)" << std::endl;
   for(unsigned int i=0; i<A.NumRows(); i++)
   {
      for(unsigned int j=0; j<A.NumCols(); j++)
         std::cout << A[i][j] << " ";
      std::cout << std::endl;
   }

   // the solution                                                                
   std::cout << "Lanczos ... " << std::endl;
   NTL::vec_ZZ_p X;
   if(!Lanczos(A,X,y))
      std::cout << "Lanczos ERROR!" << std::endl;  


   //std::cout << "Reduce matrix " << nrows << "x" << ncols << " ..." << std::endl;
   //reduce_matrix(extra_rels, &nrows, &ncols, M);
 
   //std::cout << "Block lanczos ..." << std::endl;
   //int cnt=0;
   //do
   //{                                                                            
   //   nullrows = block_lanczos(nrows, 0, ncols, M);        
   //} while (nullrows == NULL /*&& cnt++<1000*/);
   

   /*
   puts("\nNULL SPACES");
   for(uint32_t i=0; i<nrows; i++)
   {
      for(uint32_t j=0; j<ncols; j++)
      {
         if(get_null_entry(nullrows, i, j))
            printf("1 ");
         else
            printf("0 ");
      }
      puts("");
   }

   uint64_t mask;                                                                

   for (i = 0, mask = 0; i < (long)ncols; i++) // create mask of nullspace vectors 
      mask |= nullrows[i];                                                     

   int count;                                                                                
   for (i = count = 0; i < 64; i++) / count nullspace vectors found 
   {                                                                            
      if (mask & ((uint64_t)(1) << i))                                          
         count++;                                                               
   }   

   for (count = 0; count < 64; count++)                                         
   {                                                                            
      if (mask & ((uint64_t)(1) << count))                                     
      { 
         //printf("count: %d\n", count);                                                                       
      }                                                                        
   }     

   NTL::mat_GF2 tmp1 = matrix.sM;
   NTL::gauss(tmp1);
   std::cout << "NTL gauss: " << std::endl << tmp1 << std::endl;

   NTL::mat_GF2 tmp2 = matrix.sM;
   NTL::kernel(tmp2, tmp2);
   std::cout << "NTL kernel: " << std::endl << tmp2 << std::endl;

   NTL::mat_GF2 tmp3 = matrix.sM;
   NTL::image(tmp3, tmp3);
   std::cout << "NTL image: " << std::endl << tmp3 << std::endl;
   */

#endif

   
   bool stop;

   // Solve Mx = 0
   NTL::GF2 temp;

   for(i=0; i < matrix.col-1; i++)
   {
      matrix.sfreeCols[i] = 0;
      stop = false;
      j = i;
      if(matrix.sM[i][j]==0)
      {
         // find k so that M[k,j]==1
         k = i;
         while(++k < matrix.col-1 && matrix.sM[k][j]==0)
         {
         }
         // swap row i with row k
         if(k >= matrix.col-1)
         {
            matrix.sfreeCols[i] = 1;
            stop = true;
         }
         for(j = 0; j <= matrix.row-1 && !stop; j++)
         {
            temp = matrix.sM[i][j];
            matrix.sM[i][j] = matrix.sM[k][j];
            matrix.sM[k][j] = temp;
         }

      }
      for(k = 0; k < matrix.col-1 && !stop; k++)
      {
         if(matrix.sM[k][i]==1 && k != i)
         {
            for(j = 0; j <= matrix.row-1; j++)
            {
               matrix.sM[k][j] =  
               ((NTL::IsOne(matrix.sM[k][j]) ||NTL::IsOne(matrix.sM[i][j])) &&
               ((NTL::IsZero(matrix.sM[k][j])||NTL::IsZero(matrix.sM[i][j]))));
            }
         }
      }
   }
   

#ifdef LANCZOS      
   std::cout << "OLD MATRIX RESULTS" << std::endl; 
   std::cout << matrix.sM << std::endl;
   std::cout << matrix.sfreeCols << std::endl;
 
   //NTL::kernel(matrix.sM, matrix.sM);
   //NTL::gauss(matrix.sM);
   //std::cout << matrix.sM << std::endl;  
   /*
   for(uint32_t i=0; i<nrows; i++)
   {
      int cnt_ones=0;
      for(uint32_t j=0; j<ncols; j++)
      {
         if(get_null_entry(nullrows, i, j))
         {
            matrix.sfreeCols[j]=1;
            cnt_ones++;
         }
         else
         {
            matrix.sfreeCols[j]=0;
         }
      }
      
      if(cnt_ones>0) break;
   }


   std::cout << matrix.sfreeCols << std::endl;
   */
#endif

}
// }}}



