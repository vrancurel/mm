
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


/*
  File with basic error-handling operations
  Based on previous version on Zilog
  System 8000 setret() etc.
  Ported to Pyramid 9810 late 1987
  */

static	char	rcsid[] = "$Id: err.c,v 1.6 1995/01/30 14:49:14 des Exp $";

#include	<stdio.h>
#include	<setjmp.h>
#include	<ctype.h>
#include        "err.h"


/* set_err_flag -- sets err_flag -- returns old err_flag */
int	set_err_flag(flag)
int	flag;
{
   return 0;
}

/* ev_err -- reports error (err_num) in file "file" at line "line_num" and
   returns to user error handler;
   list_num is an error list number (0 is the basic list 
   pointed by err_mesg, 1 is the basic list of warnings)
 */
int	ev_err(file,err_num,line_num,fn_name,list_num)
char	*file, *fn_name;
int	err_num, line_num,list_num;
{
   return 0;
}
