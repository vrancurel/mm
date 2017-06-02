
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


/* matop.c 1.3 11/25/87 */


#include	<stdio.h>
#include	"matrix.h"

static	char	rcsid[] = "$Id: matop.c,v 1.3 1994/01/13 05:30:28 des Exp $";


/* vm_mlt -- vector-matrix multiplication 
		-- Note: b is treated as a row vector */
VEC	*vm_mlt(A,b,out)
MAT	*A;
VEC	*b,*out;
{
	u_int	j,m,n;
	/* Real	sum,**A_v,*b_v; */

	if ( A==(MAT *)NULL || b==(VEC *)NULL )
		error(E_NULL,"vm_mlt");
	if ( A->m != b->dim )
		error(E_SIZES,"vm_mlt");
	if ( b == out )
		error(E_INSITU,"vm_mlt");
	if ( out == (VEC *)NULL || out->dim != A->n )
		out = v_resize(out,A->n);

	m = A->m;		n = A->n;

	v_zero(out);
	for ( j = 0; j < m; j++ )
		if ( b->ve[j] != 0.0 )
		    __mltadd__(out->ve,A->me[j],b->ve[j],(int)n);
	/**************************************************
	A_v = A->me;		b_v = b->ve;
	for ( j=0; j<n; j++ )
	{
		sum = 0.0;
		for ( i=0; i<m; i++ )
			sum += b_v[i]*A_v[i][j];
		out->ve[j] = sum;
	}
	**************************************************/

	return out;
}
