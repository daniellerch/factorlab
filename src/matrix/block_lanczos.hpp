
#ifndef __GNFS_BLOCK_LANCZOS_HPP__
#define __GNFS_BLOCK_LANCZOS_HPP__

#include <cstdint>

#define MAX(x, y) ((x) > (y) ? (x) : (y))

typedef struct la_col_t /* matrix column */                                      
{                                                                                
   uint32_t * data;     /* The list of occupied rows in this column */              
   uint32_t weight;     /* Number of nonzero entries in this column */              
   uint32_t orig;         /* Original relation number */                            
} la_col_t; 


/*--------------------------------------------------------------------*/
static __inline__ void insert_col_entry(la_col_t * col, int32_t entry)             
{                                                                                
   if (((col->weight >> 4) << 4) == col->weight) /* need more space */           
   {                                                                             
       if (col->weight != 0) col->data =                                         
           (uint32_t *) realloc(col->data, (col->weight + 16)*sizeof(int32_t)); 
       else col->data = (uint32_t *) malloc(16*sizeof(int32_t));                
   }                                                                             
                                                                                 
   col->data[col->weight] = entry;                                               
   col->weight++;                                                                
}   

/*--------------------------------------------------------------------*/
static __inline__ void copy_col(la_col_t * col2, la_col_t * col1)                
{                                                                                
   col2->weight = col1->weight;                                                  
   col2->data = col1->data;                                                      
   col2->orig = col1->orig;                                                      
}                                                                                
                                                                                 
/*--------------------------------------------------------------------*/
static __inline__ void swap_cols(la_col_t * col2, la_col_t * col1)               
{                                                                                
   la_col_t temp;                                                                
                                                                                 
   temp.weight = col1->weight;                                                   
   temp.data = col1->data;                                                       
   temp.orig = col1->orig;                                                       
                                                                                 
   col1->weight = col2->weight;                                                  
   col1->data = col2->data;                                                      
   col1->orig = col2->orig;                                                      
                                                                                 
   col2->weight = temp.weight;                                                   
   col2->data = temp.data;                                                       
   col2->orig = temp.orig;                                                       
}    

/*--------------------------------------------------------------------*/
static __inline__ void clear_col(la_col_t * col)                                 
{                                                                                
   col->weight = 0;                                                              
}                                                                                
                                                                                 
/*--------------------------------------------------------------------*/
static __inline__ void free_col(la_col_t * col)                                  
{                                                                                
   if (col->weight) free(col->data);                                       
}  






la_col_t *create_matrix(int64_t size);
void reduce_matrix(uint32_t extra_rels, uint32_t *nrows, uint32_t *ncols, la_col_t *cols);
uint64_t * block_lanczos(uint32_t nrows, uint32_t dense_rows, uint32_t ncols, la_col_t *B);
uint64_t get_null_entry(uint64_t * nullrows, uint32_t i, uint32_t l);

#endif

