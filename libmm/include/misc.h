/* misc.h */

#ifndef __MISC_H__
#define __MISC_H__	1

/* PROTO misc.c */
/* src/misc.c */
void xperror(char *str);
void *xrealloc(void *ptr, size_t size);
void *xmalloc(size_t size);
char *xstrdup(char *str);
#endif
