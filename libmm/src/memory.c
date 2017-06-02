
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/* memory.c 1.3 11/25/87 */

#include 	"matrix.h"


static	char	rcsid[] = "$Id: memory.c,v 1.13 1994/04/05 02:10:37 des Exp $";

/* m_get -- gets an mxn matrix (in MAT form) by dynamic memory allocation */
MAT	*m_get(m,n)
int	m,n;
{
   MAT	*matrix;
   int	i;
   
   if (m < 0 || n < 0)
     error(E_NEG,"m_get");

   if ((matrix=NEW(MAT)) == (MAT *)NULL )
     error(E_MEM,"m_get");
   
   matrix->m = m;		matrix->n = matrix->max_n = n;
   matrix->max_m = m;	matrix->max_size = m*n;
#ifndef SEGMENTED
   if ((matrix->base = NEW_A(m*n,Real)) == (Real *)NULL )
   {
      free(matrix);
      error(E_MEM,"m_get");
   }
#else
   matrix->base = (Real *)NULL;
#endif
   if ((matrix->me = (Real **)calloc(m,sizeof(Real *))) == 
       (Real **)NULL )
   {	free(matrix->base);	free(matrix);
	error(E_MEM,"m_get");
     }
   
#ifndef SEGMENTED
   /* set up pointers */
   for ( i=0; i<m; i++ )
     matrix->me[i] = &(matrix->base[i*n]);
#else
   for ( i = 0; i < m; i++ )
     if ( (matrix->me[i]=NEW_A(n,Real)) == (Real *)NULL )
       error(E_MEM,"m_get");
#endif
   
   return (matrix);
}


/* px_get -- gets a PERM of given 'size' by dynamic memory allocation
   -- Note: initialized to the identity permutation */
PERM	*px_get(size)
int	size;
{
   PERM	*permute;
   int	i;

   if (size < 0)
     error(E_NEG,"px_get");

   if ((permute=NEW(PERM)) == (PERM *)NULL )
     error(E_MEM,"px_get");
   
   permute->size = permute->max_size = size;
   if ((permute->pe = NEW_A(size,u_int)) == (u_int *)NULL )
     error(E_MEM,"px_get");
   
   for ( i=0; i<size; i++ )
     permute->pe[i] = i;
   
   return (permute);
}

/* v_get -- gets a VEC of dimension 'dim'
   -- Note: initialized to zero */
VEC	*v_get(size)
int	size;
{
   VEC	*vector;
   
   if (size < 0)
     error(E_NEG,"v_get");

   if ((vector=NEW(VEC)) == (VEC *)NULL )
     error(E_MEM,"v_get");
   
   vector->dim = vector->max_dim = size;
   if ((vector->ve=NEW_A(size,Real)) == (Real *)NULL )
   {
      free(vector);
      error(E_MEM,"v_get");
   }
   
   return (vector);
}

/* m_free -- returns MAT & asoociated memory back to memory heap */
int	m_free(mat)
MAT	*mat;
{
#ifdef SEGMENTED
   int	i;
#endif
   
   if ( mat==(MAT *)NULL || (int)(mat->m) < 0 ||
       (int)(mat->n) < 0 )
     /* don't trust it */
     return (-1);
   
#ifndef SEGMENTED
   if ( mat->base != (Real *)NULL ) {
      free((char *)(mat->base));
   }
#else
   for ( i = 0; i < mat->max_m; i++ )
     if ( mat->me[i] != (Real *)NULL ) {
	free((char *)(mat->me[i]));
     }
#endif
   if ( mat->me != (Real **)NULL ) {
      free((char *)(mat->me));
   }
   
   free((char *)mat);
   
   return (0);
}



/* px_free -- returns PERM & asoociated memory back to memory heap */
int	px_free(px)
PERM	*px;
{
   if ( px==(PERM *)NULL || (int)(px->size) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( px->pe == (u_int *)NULL ) {
      free((char *)px);
   }
   else
   {
      free((char *)px->pe);
      free((char *)px);
   }
   
   return (0);
}



/* v_free -- returns VEC & asoociated memory back to memory heap */
int	v_free(vec)
VEC	*vec;
{
   if ( vec==(VEC *)NULL || (int)(vec->dim) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( vec->ve == (Real *)NULL ) {
      free((char *)vec);
   }
   else
   {
      free((char *)vec->ve);
      free((char *)vec);
   }
   
   return (0);
}



/* m_resize -- returns the matrix A of size new_m x new_n; A is zeroed
   -- if A == NULL on entry then the effect is equivalent to m_get() */
MAT	*m_resize(A,new_m,new_n)
MAT	*A;
int	new_m, new_n;
{
   int	i;
   int	new_max_m, new_max_n, new_size, old_m, old_n;
   
   if (new_m < 0 || new_n < 0)
     error(E_NEG,"m_resize");

   if ( ! A )
     return m_get(new_m,new_n);

   /* nothing was changed */
   if (new_m == A->m && new_n == A->n)
     return A;

   old_m = A->m;	old_n = A->n;
   if ( new_m > A->max_m )
   {	/* re-allocate A->me */

      A->me = RENEW(A->me,new_m,Real *);
      if ( ! A->me )
	error(E_MEM,"m_resize");
   }
   new_max_m = max(new_m,A->max_m);
   new_max_n = max(new_n,A->max_n);
   
#ifndef SEGMENTED
   new_size = new_max_m*new_max_n;
   if ( new_size > A->max_size )
   {	/* re-allocate A->base */

      A->base = RENEW(A->base,new_size,Real);
      if ( ! A->base )
	error(E_MEM,"m_resize");
      A->max_size = new_size;
   }
   
   /* now set up A->me[i] */
   for ( i = 0; i < new_m; i++ )
     A->me[i] = &(A->base[i*new_n]);
   
   /* now shift data in matrix */
   if ( old_n > new_n )
   {
      for ( i = 1; i < min(old_m,new_m); i++ )
	MEM_COPY((char *)&(A->base[i*old_n]),
		 (char *)&(A->base[i*new_n]),
		 sizeof(Real)*new_n);
   }
   else if ( old_n < new_n )
   {
      for ( i = (int)(min(old_m,new_m))-1; i > 0; i-- )
      {   /* copy & then zero extra space */
	 MEM_COPY((char *)&(A->base[i*old_n]),
		  (char *)&(A->base[i*new_n]),
		  sizeof(Real)*old_n);
	 __zero__(&(A->base[i*new_n+old_n]),(new_n-old_n));
      }
      __zero__(&(A->base[old_n]),(new_n-old_n));
      A->max_n = new_n;
   }
   /* zero out the new rows.. */
   for ( i = old_m; i < new_m; i++ )
     __zero__(&(A->base[i*new_n]),new_n);
#else
   if ( A->max_n < new_n )
   {
      Real	*tmp;
      
      for ( i = 0; i < A->max_m; i++ )
      {
	 if ( (tmp = RENEW(A->me[i],new_max_n,Real)) == NULL )
	   error(E_MEM,"m_resize");
	 else {	
	    A->me[i] = tmp;
	 }
      }
      for ( i = A->max_m; i < new_max_m; i++ )
      {
	 if ( (tmp = NEW_A(new_max_n,Real)) == NULL )
	   error(E_MEM,"m_resize");
	 else {
	    A->me[i] = tmp;

	 }
      }
   }
   else if ( A->max_m < new_m )
   {
      for ( i = A->max_m; i < new_m; i++ ) 
	if ( (A->me[i] = NEW_A(new_max_n,Real)) == NULL )
	  error(E_MEM,"m_resize");
      
   }
   
   if ( old_n < new_n )
   {
      for ( i = 0; i < old_m; i++ )
	__zero__(&(A->me[i][old_n]),new_n-old_n);
   }
   
   /* zero out the new rows.. */
   for ( i = old_m; i < new_m; i++ )
     __zero__(A->me[i],new_n);
#endif
   
   A->max_m = new_max_m;
   A->max_n = new_max_n;
   A->max_size = A->max_m*A->max_n;
   A->m = new_m;	A->n = new_n;
   
   return A;
}

/* px_resize -- returns the permutation px with size new_size
   -- px is set to the identity permutation */
PERM	*px_resize(px,new_size)
PERM	*px;
int	new_size;
{
   int	i;
   
   if (new_size < 0)
     error(E_NEG,"px_resize");

   if ( ! px )
     return px_get(new_size);
   
   /* nothing is changed */
   if (new_size == px->size)
     return px;

   if ( new_size > px->max_size )
   {
      px->pe = RENEW(px->pe,new_size,u_int);
      if ( ! px->pe )
	error(E_MEM,"px_resize");
      px->max_size = new_size;
   }
   if ( px->size <= new_size )
     /* extend permutation */
     for ( i = px->size; i < new_size; i++ )
       px->pe[i] = i;
   else
     for ( i = 0; i < new_size; i++ )
       px->pe[i] = i;
   
   px->size = new_size;
   
   return px;
}

/* v_resize -- returns the vector x with dim new_dim
   -- x is set to the zero vector */
VEC	*v_resize(x,new_dim)
VEC	*x;
int	new_dim;
{
   
   if (new_dim < 0)
     error(E_NEG,"v_resize");

   if ( ! x )
     return v_get(new_dim);

   /* nothing is changed */
   if (new_dim == x->dim)
     return x;

   if ( x->max_dim == 0 )	/* assume that it's from sub_vec */
     return v_get(new_dim);
   
   if ( new_dim > x->max_dim )
   {

      x->ve = RENEW(x->ve,new_dim,Real);
      if ( ! x->ve )
	error(E_MEM,"v_resize");
      x->max_dim = new_dim;
   }
   
   if ( new_dim > x->dim )
     __zero__(&(x->ve[x->dim]),new_dim - x->dim);
   x->dim = new_dim;
   
   return x;
}
