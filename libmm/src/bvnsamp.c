/* bvnsamp.c */

#include "mm.h"
#include <math.h>

/** 
 * sample bivariate correlated random normal variables
 * 
 * @param mu1 mean 1
 * @param sigma1 deviation 1
 * @param mu2 mean 2
 * @param sigma2 deviation 2
 * @param rho correlation factor
 * @param x1p x1 result pointer (or NULL)
 * @param x2p x2 result pointer (or NULL)
 */
void bvnsamp(double mu1, double sigma1,
	     double mu2, double sigma2, double rho,
	     double *x1p, double *x2p)
{
  double z1, z2, x1, x2;
  
  z1 = invnorm(unirandp());
  z2 = invnorm(unirandp());
  
  x1 = mu1 + sigma1*z1;
  x2 = mu2 + sigma2*(z1*rho+z2*sqrt(1-rho*rho));

  if (x1p)
    *x1p = x1;

  if (x2p)
    *x2p = x2;
}
