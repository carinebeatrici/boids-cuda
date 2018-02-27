#include <stdlib.h>
#include <stdio.h>

float **Free_matrix_real (int m_lines, int n_rows, float **matrix)
{
   int  i;  // auxiliar 
   if (matrix == NULL)return (NULL);
   if (m_lines < 1 || n_rows < 1)
     {
        printf ("** Error: Matrix too small **\n");
        return NULL;
     }
   for (i=0; i<m_lines; i++) free (matrix[i]);// free matrix lines 
   free (matrix);      // free matrix pointer 
   return (matrix);
}


float *Free_vector_real (int vector_size, float *vector)
{
   if (vector_size < 1)
     {
        printf ("** Error: vector too small **\n");
     }
   free (vector);      // free vector pointer 
   return (vector);
}



