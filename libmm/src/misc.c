/* misc.c */

#include "mm.h"

void	xperror(char *str)
{
  perror(str);
  exit(1);
}

void	*xrealloc(void *ptr, size_t size)
{
  char *nptr;

  if ((nptr = realloc(ptr, size)) == NULL)
    xperror("realloc");

  return nptr;
}

void	*xmalloc(size_t size)
{
  char *ptr;
  
  if ((ptr = malloc(size)) == NULL)
    xperror("malloc");

  return ptr;
}

char	*xstrdup(char *str)
{
  char	*n;

  if ((n = strdup(str)) == NULL)
    xperror("strdup");

  return n;
}


