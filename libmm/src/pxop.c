
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


/* pxop.c 1.5 12/03/87 */


#include	<stdio.h>
#include	"matrix.h"

static	char	rcsid[] = "$Id: pxop.c,v 1.5 1994/03/23 23:58:50 des Exp $";

/**********************************************************************
Note: A permutation is often interpreted as a matrix
		(i.e. a permutation matrix).
	A permutation px represents a permutation matrix P where
		P[i][j] == 1 if and only if px->pe[i] == j
**********************************************************************/


/* px_inv -- invert permutation -- in situ
	-- taken from ACM Collected Algorithms #250 */
PERM	*px_inv(px,out)
PERM	*px, *out;
{
    int	i, j, k, n, *p;
    
    out = px_copy(px, out);
    n = out->size;
    p = (int *)(out->pe);
    for ( n--; n>=0; n-- )
    {
	i = p[n];
	if ( i < 0 )	p[n] = -1 - i;
	else if ( i != n )
	{
	    k = n;
	    while (TRUE)
	    {
		if ( i < 0 || i >= out->size )
		    error(E_BOUNDS,"px_inv");
		j = p[i];	p[i] = -1 - k;
		if ( j == n )
		{	p[n] = i;	break;		}
		k = i;		i = j;
	    }
	}
    }
    return out;
}

/* px_vec -- permute vector */
VEC	*px_vec(px,vector,out)
PERM	*px;
VEC	*vector,*out;
{
    u_int	old_i, i, size, start;
    Real	tmp;
    
    if ( px==(PERM *)NULL || vector==(VEC *)NULL )
	error(E_NULL,"px_vec");
    if ( px->size > vector->dim )
	error(E_SIZES,"px_vec");
    if ( out==(VEC *)NULL || out->dim < vector->dim )
	out = v_resize(out,vector->dim);
    
    size = px->size;
    if ( size == 0 )
	return v_copy(vector,out);
    if ( out != vector )
    {
	for ( i=0; i<size; i++ )
	    if ( px->pe[i] >= size )
		error(E_BOUNDS,"px_vec");
	    else
		out->ve[i] = vector->ve[px->pe[i]];
    }
    else
    {	/* in situ algorithm */
	start = 0;
	while ( start < size )
	{
	    old_i = start;
	    i = px->pe[old_i];
	    if ( i >= size )
	    {
		start++;
		continue;
	    }
	    tmp = vector->ve[start];
	    while ( TRUE )
	    {
		vector->ve[old_i] = vector->ve[i];
		px->pe[old_i] = i+size;
		old_i = i;
		i = px->pe[old_i];
		if ( i >= size )
		    break;
		if ( i == start )
		{
		    vector->ve[old_i] = tmp;
		    px->pe[old_i] = i+size;
		    break;
		}
	    }
	    start++;
	}

	for ( i = 0; i < size; i++ )
	    if ( px->pe[i] < size )
		error(E_BOUNDS,"px_vec");
	    else
		px->pe[i] = px->pe[i]-size;
    }
    
    return out;
}




/* px_transp -- transpose elements of permutation
		-- Really multiplying a permutation by a transposition */
PERM	*px_transp(px,i1,i2)
PERM	*px;		/* permutation to transpose */
u_int	i1,i2;		/* elements to transpose */
{
	u_int	temp;

	if ( px==(PERM *)NULL )
		error(E_NULL,"px_transp");

	if ( i1 < px->size && i2 < px->size )
	{
		temp = px->pe[i1];
		px->pe[i1] = px->pe[i2];
		px->pe[i2] = temp;
	}

	return px;
}

/* myqsort -- a cheap implementation of Quicksort on integers
		-- returns number of swaps */
static int myqsort(a,num)
int	*a, num;
{
	int	i, j, tmp, v;
	int	numswaps;

	numswaps = 0;
	if ( num <= 1 )
		return 0;

	i = 0;	j = num;	v = a[0];
	for ( ; ; )
	{
		while ( a[++i] < v )
			;
		while ( a[--j] > v )
			;
		if ( i >= j )	break;

		tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
		numswaps++;
	}

	tmp = a[0];
	a[0] = a[j];
	a[j] = tmp;
	if ( j != 0 )
		numswaps++;

	numswaps += myqsort(&a[0],j);
	numswaps += myqsort(&a[j+1],num-(j+1));

	return numswaps;
}


/* px_sign -- compute the ``sign'' of a permutation = +/-1 where
		px is the product of an even/odd # transpositions */
int	px_sign(px)
PERM	*px;
{
	int	numtransp;
	PERM	*px2;

	if ( px==(PERM *)NULL )
		error(E_NULL,"px_sign");
	px2 = px_copy(px,PNULL);
	numtransp = myqsort(px2->pe,px2->size);
	px_free(px2);

	return ( numtransp % 2 ) ? -1 : 1;
}

