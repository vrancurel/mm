/* mvnsamp.c */

#include "mm.h"

/** 
 * generate n_dim correlated numbers using Cholesky decomposition.
 * check rho 
 * 
 * @param n_dim number of dimensions
 * @param mu mean vector
 * @param sigma deviation vector
 * @param rho correlation matrix (all values must be provided)
 * @param x output vector
 *
 * @return 0 if OK, else -1
 */
int mvnsamp(u_int n_dim, VEC *mu, VEC *sigma, MAT *rho, VEC *x)
{
  int ret;
  u_int i;
  VEC *z;
  MAT *cov;

  ret = rho_check(rho);
  assert(-1 != ret);

  cov = m_get(n_dim, n_dim);
  assert(NULL != cov);

  comp_cov(n_dim, sigma, rho, cov);

  chol_decomp(cov);

  z = v_get(n_dim);
  assert(NULL != z);

  for (i = 0;i < n_dim;i++)
    z->ve[i] = invnorm(unirandp());

  vm_mlt(cov, z, x);

  for (i = 0;i < n_dim;i++)
    x->ve[i] += mu->ve[i];

  V_FREE(z);
  M_FREE(cov);

  return 0;
}
