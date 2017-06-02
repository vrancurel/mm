
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


/* vecop.c 1.3 8/18/87 */

#include	<stdio.h>
#include	"matrix.h"

static	char	rcsid[] = "$Id: vecop.c,v 1.4 1994/03/08 05:50:39 des Exp $";

/* v_sub -- vector subtraction -- may be in-situ */
VEC	*v_sub(vec1,vec2,out)
VEC	*vec1,*vec2,*out;
{
	/* u_int	i, dim; */
	/* Real	*out_ve, *vec1_ve, *vec2_ve; */

	if ( vec1==(VEC *)NULL || vec2==(VEC *)NULL )
		error(E_NULL,"v_sub");
	if ( vec1->dim != vec2->dim )
		error(E_SIZES,"v_sub");
	if ( out==(VEC *)NULL || out->dim != vec1->dim )
		out = v_resize(out,vec1->dim);

	__sub__(vec1->ve,vec2->ve,out->ve,(int)(vec1->dim));
	/************************************************************
	dim = vec1->dim;
	out_ve = out->ve;	vec1_ve = vec1->ve;	vec2_ve = vec2->ve;
	for ( i=0; i<dim; i++ )
		out->ve[i] = vec1->ve[i]-vec2->ve[i];
		(*out_ve++) = (*vec1_ve++) - (*vec2_ve++);
	************************************************************/

	return (out);
}
