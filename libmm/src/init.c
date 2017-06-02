
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


/*
	This is a file of routines for zero-ing, and initialising
	vectors, matrices and permutations.
	This is to be included in the matrix.a library
*/

static	char	rcsid[] = "$Id: init.c,v 1.6 1994/01/13 05:36:58 des Exp $";

#include	<stdio.h>
#include	"matrix.h"

/* v_zero -- zero the vector x */
VEC	*v_zero(x)
VEC	*x;
{
	if ( x == VNULL )
		error(E_NULL,"v_zero");

	__zero__(x->ve,x->dim);
	/* for ( i = 0; i < x->dim; i++ )
		x->ve[i] = 0.0; */

	return x;
}


/* m_zero -- zero the matrix A */
MAT	*m_zero(A)
MAT	*A;
{
	int	i, A_m, A_n;
	Real	**A_me;

	if ( A == MNULL )
		error(E_NULL,"m_zero");

	A_m = A->m;	A_n = A->n;	A_me = A->me;
	for ( i = 0; i < A_m; i++ )
		__zero__(A_me[i],A_n);
		/* for ( j = 0; j < A_n; j++ )
			A_me[i][j] = 0.0; */

	return A;
}
