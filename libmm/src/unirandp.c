/* unirandp.c */

#include "mm.h"
#include <math.h>

/** 
 * @return a uniform random probability number between 0.0 and 1.0
 */
double unirandp(void)
{
  double n;
  
  n = (double)rand()/(double)RAND_MAX;

  return n;
}
