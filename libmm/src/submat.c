
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


/* 1.2 submat.c 11/25/87 */

#include	<stdio.h>
#include	"matrix.h"

static	char	rcsid[] = "$Id: submat.c,v 1.2 1994/01/13 05:28:12 des Exp $";

/* _set_col -- sets column of matrix to values given in vec (in situ) */
MAT	*_set_col(mat,col,vec,i0)
MAT	*mat;
VEC	*vec;
u_int	col,i0;
{
   u_int	i,lim;
   
   if ( mat==(MAT *)NULL || vec==(VEC *)NULL )
     error(E_NULL,"_set_col");
   if ( col >= mat->n )
     error(E_RANGE,"_set_col");
   lim = min(mat->m,vec->dim);
   for ( i=i0; i<lim; i++ )
     mat->me[i][col] = vec->ve[i];
   
   return (mat);
}
