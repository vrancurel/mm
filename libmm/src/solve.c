
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
	Matrix factorisation routines to work with the other matrix files.
*/

/* solve.c 1.2 11/25/87 */
static	char	rcsid[] = "$Id: solve.c,v 1.3 1994/01/13 05:29:57 des Exp $";

#include	<stdio.h>
#include        "matrix2.h"
#include	<math.h>





/* Most matrix factorisation routines are in-situ unless otherwise specified */

/* Usolve -- back substitution with optional over-riding diagonal
		-- can be in-situ but doesn't need to be */
VEC	*Usolve(matrix,b,out,diag)
MAT	*matrix;
VEC	*b, *out;
double	diag;
{
	u_int	dim /* , j */;
	int	i, i_lim;
	Real	**mat_ent, *mat_row, *b_ent, *out_ent, *out_col, sum, tiny;

	if ( matrix==(MAT *)NULL || b==(VEC *)NULL )
		error(E_NULL,"Usolve");
	dim = min(matrix->m,matrix->n);
	if ( b->dim < dim )
		error(E_SIZES,"Usolve");
	if ( out==(VEC *)NULL || out->dim < dim )
		out = v_resize(out,matrix->n);
	mat_ent = matrix->me;	b_ent = b->ve;	out_ent = out->ve;

	tiny = 10.0/HUGE_VAL;

	for ( i=dim-1; i>=0; i-- )
		if ( b_ent[i] != 0.0 )
		    break;
		else
		    out_ent[i] = 0.0;
	i_lim = i;

	for (    ; i>=0; i-- )
	{
		sum = b_ent[i];
		mat_row = &(mat_ent[i][i+1]);
		out_col = &(out_ent[i+1]);
		sum -= __ip__(mat_row,out_col,i_lim-i);
		/******************************************************
		for ( j=i+1; j<=i_lim; j++ )
			sum -= mat_ent[i][j]*out_ent[j];
			sum -= (*mat_row++)*(*out_col++);
		******************************************************/
		if ( diag==0.0 )
		{
			if ( fabs(mat_ent[i][i]) <= tiny*fabs(sum) )
				error(E_SING,"Usolve");
			else
				out_ent[i] = sum/mat_ent[i][i];
		}
		else
			out_ent[i] = sum/diag;
	}

	return (out);
}

/* Lsolve -- forward elimination with (optional) default diagonal value */
VEC	*Lsolve(matrix,b,out,diag)
MAT	*matrix;
VEC	*b,*out;
double	diag;
{
	u_int	dim, i, i_lim /* , j */;
	Real	**mat_ent, *mat_row, *b_ent, *out_ent, *out_col, sum, tiny;

	if ( matrix==(MAT *)NULL || b==(VEC *)NULL )
		error(E_NULL,"Lsolve");
	dim = min(matrix->m,matrix->n);
	if ( b->dim < dim )
		error(E_SIZES,"Lsolve");
	if ( out==(VEC *)NULL || out->dim < dim )
		out = v_resize(out,matrix->n);
	mat_ent = matrix->me;	b_ent = b->ve;	out_ent = out->ve;

	for ( i=0; i<dim; i++ )
		if ( b_ent[i] != 0.0 )
		    break;
		else
		    out_ent[i] = 0.0;
	i_lim = i;

	tiny = 10.0/HUGE_VAL;

	for (    ; i<dim; i++ )
	{
		sum = b_ent[i];
		mat_row = &(mat_ent[i][i_lim]);
		out_col = &(out_ent[i_lim]);
		sum -= __ip__(mat_row,out_col,(int)(i-i_lim));
		/*****************************************************
		for ( j=i_lim; j<i; j++ )
			sum -= mat_ent[i][j]*out_ent[j];
			sum -= (*mat_row++)*(*out_col++);
		******************************************************/
		if ( diag==0.0 )
		{
			if ( fabs(mat_ent[i][i]) <= tiny*fabs(sum) )
				error(E_SING,"Lsolve");
			else
				out_ent[i] = sum/mat_ent[i][i];
		}
		else
			out_ent[i] = sum/diag;
	}

	return (out);
}

