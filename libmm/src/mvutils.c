/* mvutils.c */

#include "mm.h"

void mat_print(MAT *mat)
{
  u_int i, j;

  for (i = 0;i < mat->m;i++)
    {
      printf("|");
      for (j = 0;j < mat->n;j++)
        printf("%16f ", mat->me[i][j]);
      printf("|\n");
    }
}

void perm_print(PERM *perm)
{
  u_int i, j;

  for (i = 0;i < perm->size;i++)
    {
      printf("|");
      for (j = 0;j < perm->size;j++)
        {
          if (j == perm->pe[i])
            printf("1 ");
          else
            printf("0 ");
        }
      printf("|\n");
    }
}

void vec_print(VEC *vec)
{
  int i;

  for (i = 0;i < vec->dim;i++)
    printf("|%16f|\n", vec->ve[i]);
}

/*
 * I/O
 */

/** 
 * use scilab format
 * 
 * @param mat 
 */
void mat_print_scilab(MAT *mat, FILE *f)
{
  u_int i, j;

  fprintf(f, "[");
  for (i = 0;i < mat->m;i++)
    {
      for (j = 0;j < mat->n;j++)
        fprintf(f, "%f ", mat->me[i][j]);
      if (i != (mat->m - 1))
        fprintf(f, ";");
    }
  fprintf(f, "]\n");
}

void mat_parse_scilab(MAT *mat, FILE *f)
{
  int start = 0;
  char buf[256];
  int i = 0;
  int j = 0;
  int k = 0;

  while (1)
    {
      int c = fgetc(f);
      
      switch (c)
        {
        case '[':
          start = 1;
          break ;
        case ' ':
          if (!start)
            continue ;
          buf[k] = 0;
          mat->me[i][j] = atof(buf);
          j++;
          k = 0;
          break ;
        case ';':
          i++;
          j = 0;
          break ;
        case ']':
          return ;
        default:
          buf[k++] = c;
          break ;
        }
    }
}

void vec_print_scilab(VEC *vec, FILE *f)
{
  int i;

  fprintf(f, "[");
  for (i = 0;i < vec->dim;i++)
    fprintf(f, "%f ", vec->ve[i]);
  fprintf(f, "]\n");
}

void vec_parse_scilab(VEC *vec, FILE *f)
{
  int start = 0;
  char buf[256];
  int i = 0;
  int k = 0;

  while (1)
    {
      int c = fgetc(f);
      
      switch (c)
        {
        case '[':
          start = 1;
          break ;
        case ' ':
          if (!start)
            continue ;
          buf[k] = 0;
          vec->ve[i] = atof(buf);
          i++;
          k = 0;
          break ;
        case ']':
          return ;
        default:
          buf[k++] = c;
          break ;
        }
    }
}

/*
 *
 */

int mat_init(MAT *mat, double *dat, u_int n)
{
  u_int i, j;

  assert(n == mat->m*mat->n);

  for (i = 0;i < mat->m;i++)
    for (j = 0;j < mat->n;j++)
      mat->me[i][j] = dat[i*mat->m+j];

  return 0;
}

int vec_init(VEC *vec, double *dat, u_int n)
{
  u_int i;

  assert(n == vec->dim);

  for (i = 0;i < n;i++)
    vec->ve[i] = dat[i];

  return 0;
}

/** 
 * check correlation matrix
 * 
 * @param rho 
 * 
 * @return 
 */
int rho_check(MAT *rho)
{
  u_int i, j;

  assert(rho->m == rho->n);

  for (i = 0;i < rho->m;i++)
    {
      assert(rho->me[i][i] == 1.);
    }

  for (i = 0;i < rho->m;i++)
    for (j = 0;j < i;j++)
      {
        assert(rho->me[i][j] == rho->me[j][i]);
      }

  return 0;
}

/** 
 * compute the covariance matrix
 * 
 * @param n_dim number of dimensions
 * @param sigma deviation vector
 * @param rho correlation matrix (not checked)
 * @param cov covariance matrix
 */
void comp_cov(u_int n_dim, VEC *sigma, MAT *rho, MAT *cov)
{
  u_int i, j;

  for (i = 0;i < n_dim;i++)
    for (j = 0;j < n_dim;j++)
      cov->me[i][j] = rho->me[i][j] * sigma->ve[i] * sigma->ve[j];
}

/** 
 * compute the correlation matrix
 * 
 * @param n_dim 
 * @param cov 
 * @param rho 
 */
void comp_rho(u_int n_dim, MAT *cov, VEC *sigma, MAT *rho)
{
  u_int i, j;

  for (i = 0;i < n_dim;i++)
    for (j = 0;j < n_dim;j++)
      {
        double prod;

        prod = sigma->ve[i]*sigma->ve[j];

        assert(0 != prod);

        rho->me[i][j] = cov->me[i][j]/prod;
      }
}

/** 
 * perform a Cholesky decomposition (of the covariance matrix)
 * 
 * @param cov 
 */
void chol_decomp(MAT *cov)
{
  u_int i, j;

  CHfactor(cov);

  /*
   * reset lower triangle to zero
   */
  for (i = 1;i < cov->m;i++)
    for (j = 0; j < i;j++)
      cov->me[i][j] = 0.0;
}

/** 
 * non recursive factorial
 * 
 * @param n 
 * 
 * @return 
 */
int nr_fact(int n)
{
  int i;

  if (0 == n)
    return 1;

  for (i = n-1;i > 0;i--)
    n *= i;
  
  return n;
}

/** 
 * generate kth permutation
 * 
 * @param k 0<= k =n!
 * @param perm must contain values from 1 to n
 */
void gen_perm(int k, PERM *sigma)
{
  u_int j;
  
  for (j = 2;j < sigma->size + 1;j++)
    {
      u_int tmp;
      
      k = k/(j-1);
      tmp = sigma->pe[(k%j)];
      sigma->pe[(k%j)] = sigma->pe[j-1];
      sigma->pe[j-1] = tmp;
    }
}

/** 
 * return the determinant of a triangular (square) matrix (product of
 * diagonal elts). No checks are made on matrix
 * 
 * @param mat 
 * 
 * @return 
 */
double diag_prod(MAT *mat)
{
  u_int i;
  double prod;

  prod = 1.0;
  for (i = 0;i < mat->m;i++)
    prod *= mat->me[i][i];

  return prod;
}

/** 
 * compute determinant in using LUP decomposition
 * 
 * @param n_dim 
 * @param mat 
 * @param result 
 * 
 * @return 
 */
int comp_det(u_int n_dim, MAT *mat, double *result)
{
  MAT *m, *l, *u;
  PERM *p;
  u_int i, j;

  m = NULL;
  m = m_copy(mat, m);
  assert(NULL != m);

  p = NULL;
  p = px_resize(p, mat->m);
  assert(NULL != p);

  LUfactor(m, p);

  l = NULL;
  l = m_copy(m, l);
  assert(NULL != l);

  u = NULL;
  u = m_copy(m, u);
  assert(NULL != u);

  for (i = 1;i < n_dim;i++)
    for (j = 0;j < i;j++)
      u->me[i][j] = 0;

  for (j = 0;j < n_dim;j++)
    for (i = 0;i <= j;i++)
      l->me[i][j] = (i == j) ? 1 : 0;

  px_inv(p, p);

  /*
   * now det(mat) = det(P^-1) * det(l) * det(u)
   */
  *result = px_sign(p) * diag_prod(l) * diag_prod(u);

  M_FREE(u);
  M_FREE(l);
  PX_FREE(p);
  M_FREE(m);

  return 0;
}
