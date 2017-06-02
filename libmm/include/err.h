
/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
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


/* err.h  28/09/1993 */

/*  RCS id: $Id: err.h,v 1.2 1995/01/30 14:48:05 des Exp $  */


#ifndef ERRHEADER
#define ERRHEADER


#include        "machine.h"

/* Error recovery */

extern	int ev_err(char *,int,int,char *,int);  /* main error handler */
extern	int set_err_flag(int flag);         /* for different ways of handling */


/* error(E_TYPE,"myfunc") raises error type E_TYPE for function my_func() */
#define	error(err_num,fn_name)	ev_err(__FILE__,err_num,__LINE__,fn_name,0)

/* warning(WARN_TYPE,"myfunc") raises warning type WARN_TYPE for 
   function my_func() */
#define warning(err_num,fn_name) ev_err(__FILE__,err_num,__LINE__,fn_name,1) 


/* error types */
#define	E_UNKNOWN	0
#define	E_SIZES		1
#define	E_BOUNDS	2
#define	E_MEM		3
#define	E_SING		4
#define	E_POSDEF	5
#define	E_FORMAT	6
#define	E_INPUT		7
#define	E_NULL		8
#define	E_SQUARE	9
#define	E_RANGE		10
#define	E_INSITU2	11
#define	E_INSITU	12
#define	E_ITER		13
#define	E_CONV		14
#define	E_START		15
#define	E_SIGNAL	16
#define	E_INTERN	17
#define	E_EOF		18
#define E_SHARED_VECS   19
#define E_NEG           20
#define E_OVERWRITE     21
#define E_BREAKDOWN     22

/* warning types */
#define WARN_UNKNOWN     	0
#define WARN_WRONG_TYPE 	1
#define WARN_NO_MARK		2
#define WARN_RES_LESS_0         3
#define WARN_SHARED_VEC		4


#endif   /* ERRHEADER */

