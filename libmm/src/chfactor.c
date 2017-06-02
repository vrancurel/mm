
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

/* CHfactor.c 1.2 11/25/87 */
static	char	rcsid[] = "$Id: chfactor.c,v 1.2 1994/01/13 05:36:36 des Exp $";

#include	<stdio.h>
#include	"matrix.h"
#include        "matrix2.h"
#include	<math.h>


/* Most matrix factorisation routines are in-situ unless otherwise specified */

/* CHfactor -- Cholesky L.L' factorisation of A in-situ */
MAT	*CHfactor(A)
MAT	*A;
{
	u_int	i, j, k, n;
	Real	**A_ent, *A_piv, *A_row, sum, tmp;

	if ( A==(MAT *)NULL )
		error(E_NULL,"CHfactor");
	if ( A->m != A->n )
		error(E_SQUARE,"CHfactor");
	n = A->n;	A_ent = A->me;

	for ( k=0; k<n; k++ )
	{	
		/* do diagonal element */
		sum = A_ent[k][k];
		A_piv = A_ent[k];
		for ( j=0; j<k; j++ )
		{
			/* tmp = A_ent[k][j]; */
			tmp = *A_piv++;
			sum -= tmp*tmp;
		}
		if ( sum <= 0.0 )
			error(E_POSDEF,"CHfactor");
		A_ent[k][k] = sqrt(sum);

		/* set values of column k */
		for ( i=k+1; i<n; i++ )
		{
			sum = A_ent[i][k];
			A_piv = A_ent[k];
			A_row = A_ent[i];
			sum -= __ip__(A_row,A_piv,(int)k);
			/************************************************
			for ( j=0; j<k; j++ )
				sum -= A_ent[i][j]*A_ent[k][j];
				sum -= (*A_row++)*(*A_piv++);
			************************************************/
			A_ent[j][i] = A_ent[i][j] = sum/A_ent[k][k];
		}
	}

	return (A);
}


