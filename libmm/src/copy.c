
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


static	char	rcsid[] = "$Id: copy.c,v 1.2 1994/01/13 05:37:14 des Exp $";
#include	<stdio.h>
#include	"matrix.h"


/* _m_copy -- copies matrix into new area */
MAT	*_m_copy(in,out,i0,j0)
MAT	*in,*out;
u_int	i0,j0;
{
	u_int	i /* ,j */;

	if ( in==MNULL )
		error(E_NULL,"_m_copy");
	if ( in==out )
		return (out);
	if ( out==MNULL || out->m < in->m || out->n < in->n )
		out = m_resize(out,in->m,in->n);

	for ( i=i0; i < in->m; i++ )
		MEM_COPY(&(in->me[i][j0]),&(out->me[i][j0]),
				(in->n - j0)*sizeof(Real));
		/* for ( j=j0; j < in->n; j++ )
			out->me[i][j] = in->me[i][j]; */

	return (out);
}

/* _v_copy -- copies vector into new area */
VEC	*_v_copy(in,out,i0)
VEC	*in,*out;
u_int	i0;
{
	/* u_int	i,j; */

	if ( in==VNULL )
		error(E_NULL,"_v_copy");
	if ( in==out )
		return (out);
	if ( out==VNULL || out->dim < in->dim )
		out = v_resize(out,in->dim);

	MEM_COPY(&(in->ve[i0]),&(out->ve[i0]),(in->dim - i0)*sizeof(Real));
	/* for ( i=i0; i < in->dim; i++ )
		out->ve[i] = in->ve[i]; */

	return (out);
}

/* px_copy -- copies permutation 'in' to 'out' */
PERM	*px_copy(in,out)
PERM	*in,*out;
{
	/* int	i; */

	if ( in == PNULL )
		error(E_NULL,"px_copy");
	if ( in == out )
		return out;
	if ( out == PNULL || out->size != in->size )
		out = px_resize(out,in->size);

	MEM_COPY(in->pe,out->pe,in->size*sizeof(u_int));
	/* for ( i = 0; i < in->size; i++ )
		out->pe[i] = in->pe[i]; */

	return out;
}

