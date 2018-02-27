#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float **Alocate_matrix_real (int m_lines, int n_rows)
{
   int   i;    /* auxiliar*/
   float **matrix;
   if (m_lines<1 || n_rows<1)
     {
        printf ("** Error: Matrix must be bigger than a number **\n");
        return (NULL);
     }
   /* line allocation */
   matrix = (float **) malloc (m_lines*sizeof(float*));
   if (matrix == NULL)
     {
        printf ("** Error: Not enough memory **");
        return (NULL);
     }
   /* row allocation */
   for (i=0; i<m_lines; i++)
     {
        matrix[i] = (float*) malloc (n_rows*sizeof(float));
        if (matrix[i] == NULL)
          {
             printf ("** Erro: Memoria Insuficiente **");
             return (NULL);
          }
     }
   return (matrix);
}


float *Alocate_vector_real (int vector_size)
{
   float *vector;
   if (vector_size<1)
     {
        printf ("** Error: Vector must be bigger than a number **\n");
        return (NULL);
     }
   vector = (float *) malloc (vector_size * sizeof(float));
   if (vector == NULL)
     {
        printf ("** Error: Not enough memory **");
        return (NULL);
     }
   return(vector);
}


