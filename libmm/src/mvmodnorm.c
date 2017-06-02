/* mvmodnorm.c */

#include "mm.h"
#include <math.h>

static int initfunc(tmvmix *mix, tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);

  norm->mu = v_get(mix->n_dim);
  assert(NULL != norm->mu);

  norm->cov = m_get(mix->n_dim, mix->n_dim);
  assert(NULL != norm->cov);

  return 0;
}

static void resetfunc(tmvmix *mix, tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);

  v_zero(norm->mu);
  m_zero(norm->cov);
}

static void copyfunc(tmvmix *mix, tmvmodinst *new, tmvmodinst *old)
{
 tmvmodnormparams *newnorm = (tmvmodnormparams *)(new + 1);
 tmvmodnormparams *oldnorm = (tmvmodnormparams *)(old + 1);

 v_copy(oldnorm->mu, newnorm->mu);
 m_copy(oldnorm->cov, newnorm->cov);
}

static void randfunc(tmvmix *mix, tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);
  u_int i;

  m_zero(norm->cov);

  for (i = 0;i < mix->n_dim;i++)
    {
      double length;

      length = abs(mix->ranges[i].end - mix->ranges[i].start);
      
      norm->mu->ve[i] = mix->ranges[i].start + length * unirandp();
      norm->cov->me[i][i] = (length/4 + 1.0)*(length/4 + 1.0);
    }
}

struct normdd
{
  MAT *invcov;
  double det;
};

static int densityinitfunc(tmvmix *mix, tmvmodinst *model, void **ddp)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);
  double det;
  struct normdd *dd;
  int ret;

#if 0
  printf("covariance matrix:\n");
  mat_print(norm->cov);
#endif

  ret = comp_det(mix->n_dim, norm->cov, &det);
#if 0
  printf("determinant: %g\n", det);
#endif
  assert(-1 != ret);

  //singular matrix, inverse not computable
  //assert(0.0 != det);
  if (0.0 == det)
    return -1;

  if (det < 0.0)
    return -1;

  dd = xmalloc(sizeof (*dd));

  dd->det = det;

  dd->invcov = m_get(mix->n_dim, mix->n_dim);
  assert(NULL != dd->invcov);

  m_inverse(norm->cov, dd->invcov);

  *ddp = dd;

  return 0;
}

static void densityfreefunc(tmvmix *mix, tmvmodinst *model, void *ddp)
{
  struct normdd *dd = (struct normdd *)ddp;

  M_FREE(dd->invcov);
  free(dd);
}

static int densityfunc(tmvmix *mix,
                       tmvmodinst *model, 
                       void *ddp,
                       VEC *sample, 
                       double *result)
{
  struct normdd *dd = (struct normdd *)ddp;
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);
  double k;
  VEC *csample, *r1;
  double r2;
  u_int j, r;

  k = mix->n_dim;

  csample = v_get(mix->n_dim);
  assert(NULL != csample);

  v_sub(sample, norm->mu, csample); //centered sample

  r1 = v_get(mix->n_dim);
  assert(NULL != r1);

  /*
   * compute r1 = csample^T * invcov
   */
  v_zero(r1);
  //for (i = 0;i < 1;i++)
  for (j = 0;j < csample->dim;j++)
    for (r = 0;r < csample->dim;r++)
      r1->ve/*[i]*/[j] += csample->ve/*[i]*/[r]*dd->invcov->me[r][j];

  /*
   * compute r2 = r1 * csample
   */
  r2 = 0.0;
  for (r = 0;r < r1->dim;r++)
    r2 += r1->ve[r]*csample->ve[r];

  *result = 1./(pow(2.*M_PI,k/2.)*sqrt(dd->det))*exp((-1./2)*r2);

  V_FREE(r1);
  V_FREE(csample);

  return 0;
}

static int jointcdffunc(tmvmix *mix,
                        tmvmodinst *model, 
                        void *ddp,
                        VEC *sample, 
                        double *result)
{
  // FIXME
  return -1;
}


static int updatefunc(tmvmix *mix,
                      tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);
  u_int i, j, k;
  VEC *tmpmu;
  MAT *tmpcov;
  double tmp;

  tmpmu = v_get(mix->n_dim);
  assert(NULL != tmpmu);
  tmpcov = m_get(mix->n_dim, mix->n_dim);
  assert(NULL != tmpcov);

  v_zero(tmpmu);
  m_zero(tmpcov);

  for (i = 0;i < mix->n_samples;i++)
    for (j = 0;j < mix->n_dim;j++)
      tmpmu->ve[j] += model->probs[i] * mix->data->me[j][i];

  tmp = mix->n_samples * model->weight;

  for (i = 0;i < mix->n_dim;i++)
    norm->mu->ve[i] = tmpmu->ve[i]/tmp;

  //printf("mean found:\n");
  //vec_print(norm->mu);

  /*
   * estimation the covariance matrix
   */
  for (i = 0;i < mix->n_dim;i++)
    {
      for (j = 0;j < mix->n_dim;j++)
        {
          tmpcov->me[i][j] = 0.0;
          for (k = 0;k < mix->n_samples;k++)
            {
              tmpcov->me[i][j] += 
                model->probs[k] * 
                (mix->data->me[i][k]-norm->mu->ve[i])*
                (mix->data->me[j][k]-norm->mu->ve[j]);
            }
          norm->cov->me[i][j] = tmpcov->me[i][j]/tmp;
          if (i == j &&
              norm->cov->me[i][j] < 0.000001)
            {
              norm->cov->me[i][j] = 0.000001;
            }
        }
    }

  //printf("tmp covariance matrix found:\n");
  //mat_print(tmpcov);
  
  //printf("covariance matrix found:\n");
  //mat_print(norm->cov);

  M_FREE(tmpcov);
  V_FREE(tmpmu);

  return 0;
}

static void printfunc(tmvmix *mix,
                      tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);
  
  printf("mu:\n");
  vec_print(norm->mu);

  printf("covariance matrix:\n");
  mat_print(norm->cov);
}

static void redrawfunc(tmvmix *mix,
                       tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);

  if (mix->n_dim == 1)
    {
      //FIXME
    }
  else if (mix->n_dim == 2)
    {
      gdraw("%f+cos(t)*%f,%f+sin(t)*%f", 
            norm->mu->ve[0], sqrt(norm->cov->me[0][0]),
            norm->mu->ve[1], sqrt(norm->cov->me[1][1]));
    }
  else if (mix->n_dim == 3)
    {
      /*
       * use Mercator parametrization
       */
      //FIXME represent correlation matrix
      gdraw("%f+%f*sech(v)*cos(u),%f+%f*sech(v)*sin(u),%f+%f*tanh(v)", 
            norm->mu->ve[0], sqrt(norm->cov->me[0][0]),
            norm->mu->ve[1], sqrt(norm->cov->me[1][1]),
            norm->mu->ve[2], sqrt(norm->cov->me[2][2]));
    }
}

static void savefunc(tmvmix *mix,
                     tmvmodinst *model,
                     FILE *f)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);

  vec_print_scilab(norm->mu, f);
  mat_print_scilab(norm->cov, f);
}

static void restorefunc(tmvmix *mix,
                        tmvmodinst *model,
                        FILE *f)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);

  vec_parse_scilab(norm->mu, f);
  mat_parse_scilab(norm->cov, f);
}

static void destroyfunc(tmvmix *mix,
                        tmvmodinst *model)
{
  tmvmodnormparams *norm = (tmvmodnormparams *)(model + 1);

  V_FREE(norm->mu);
  M_FREE(norm->cov);
}

tmvmodclass mvmodnormclass = 
  {
    .name = "normal",
    .size = sizeof (tmvmodnormparams),
    .init = initfunc,
    .reset = resetfunc,
    .copy = copyfunc,
    .rand = randfunc,
    .densityinit = densityinitfunc,
    .density = densityfunc,
    .jointcdf = jointcdffunc,
    .densityfree = densityfreefunc,
    .update = updatefunc,
    .print = printfunc,
    .redraw = redrawfunc,
    .save = savefunc,
    .restore = restorefunc,
    .destroy = destroyfunc,
  };
