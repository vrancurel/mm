/* prob_test.c */

#include "mm.h"
#include <unistd.h>
#include <sys/time.h>
#include <math.h>

/*
 * general matrix utilities tests
 */

/*
 * multivariate mixture model tests
 */

#define DATFILE "/tmp/mvmix.dat"

static void test_mvmix()
{
  //#define UNIFORM_DISTR
  tmvmix *mix;
#define N_DIM 2
#define LENGTH 500
  tmvrange ranges[N_DIM];
#define N_SAMPLES	300
  

#define N_CLASSES 4
  const tmvmodclass *classes[N_CLASSES] = 
      {
        &mvmoduniclass,
        &mvmodnormclass,
        &mvmodnormclass, 
        &mvmodnormclass,
      };
  double init_weights[N_CLASSES] = 
    {
      0.1,
      1.0,
      1.0,
      1.0,
    };
  u_int i;
#ifdef UNIFORM_DISTR
  u_int j;
#endif
  FILE *f;
  MAT *data;

  gdraw("set parametric\n");
  gflush();
  
  f = fopen(DATFILE, "w+");
  assert(NULL != f);

  for (i = 0;i < N_DIM;i++)
    {
      ranges[i].start = 0;
      ranges[i].end = LENGTH;
    }

  data = m_get(N_DIM, N_SAMPLES);
  assert(NULL != data);

#ifdef UNIFORM_DISTR
  for (i = 0;i < N_SAMPLES;i++)
    {
      for (j = 0;j < N_DIM;j++)
        {
          data->me[j][i] = LENGTH*unirandp();
        }
      fprintf(f, "%f %f\n", data->me[0][i], data->me[1][i]);
    }
#else
  for (i = 0;i < N_SAMPLES/3;i++)
    {
      bvnsamp(100, 20, 100, 20, 0.34, &data->me[0][3*i], &data->me[1][3*i]);
      bvnsamp(200, 50, 300, 10, -0.21, &data->me[0][3*i+1], &data->me[1][3*i+1]);
      bvnsamp(400, 15, 400, 15, 0.15, &data->me[0][3*i+2], &data->me[1][3*i+2]);
      fprintf(f, "%f %f\n", data->me[0][3*i], data->me[1][3*i]);
      fprintf(f, "%f %f\n", data->me[0][3*i+1], data->me[1][3*i+1]);
      fprintf(f, "%f %f\n", data->me[0][3*i+2], data->me[1][3*i+2]);
    }
#endif

  (void)fclose(f);

  mix = mvmixcreate(N_DIM,
                    ranges,
                    N_SAMPLES,
                    data,
                    N_CLASSES, 
                    classes, 
                    init_weights);

  mvmixrand(mix);

  mvmixredraw(mix, DATFILE);
  getchar();

  for (i = 0;i < 100;i++)
    {
      mvmixiter(mix);
      
      mvmixredraw(mix, DATFILE);
      getchar();
    }
  
  mvmixdestroy(mix);
#undef N_DIM
#undef UNIFORM_DISTR
#undef N_SAMPLES
#undef N_CLASSES
#undef DATFILE
}

#define DATFILE "/tmp/mvmix1bis.dat"


static void test_mvmix2()
{
  //#define UNIFORM_DISTR
  tmvmix *mix;
#define N_DIM 3
#define LENGTH 500
  tmvrange ranges[N_DIM];
#define N_SAMPLES	900
  MAT *data;
#define N_CLASSES 4
  const tmvmodclass *classes[N_CLASSES] = 
      {
        &mvmoduniclass,
        &mvmodnormclass,
        &mvmodnormclass, 
        &mvmodnormclass,
      };
  double init_weights[N_CLASSES] = 
    {
      0.1,
      1.0,
      1.0,
      1.0,
    };
  u_int i;
#ifdef UNIFORM_DISTR
  u_int j;
#endif
  FILE *f;
#ifndef UNIFORM_DISTR
  int ret;
  VEC *mu1, *sigma1, *mu2, *sigma2, *mu3, *sigma3, *x;
  MAT *rho1, *rho2, *rho3;
  double mu1_dat[N_DIM] = { 100, 90, 120 };
  double sigma1_dat[N_DIM] = { 30, 20, 70 };
  double rho1_dat[N_DIM*N_DIM] =
    {
      1.0,  0.26, 0.31,
      0.26, 1.0,  0.5,
      0.31, 0.5,  1.0
    };
  double mu2_dat[N_DIM] = { 300, 50, 320 };
  double sigma2_dat[N_DIM] = { 10, 10, 25 };
  double rho2_dat[N_DIM*N_DIM] =
    {
      1.0,  0.25, 0.19,
      0.25, 1.0,  0.42,
      0.19, 0.42, 1.0
    };
  double mu3_dat[N_DIM] = { 320, 450, 350 };
  double sigma3_dat[N_DIM] = { 70, 20, 55 };
  double rho3_dat[N_DIM*N_DIM] =
    {
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0
    };
#endif
      
  gdraw("set parametric\n");
  gdraw("set isosamples 30\n");
  gdraw("sech(x) = 1/cosh(x)\n");
  gflush();
  
  f = fopen(DATFILE, "w+");
  assert(NULL != f);

  for (i = 0;i < N_DIM;i++)
    {
      ranges[i].start = 0;
      ranges[i].end = LENGTH;
    }

  data = m_get(N_DIM, N_SAMPLES);
  assert(NULL != data);

#ifdef UNIFORM_DISTR
  for (i = 0;i < N_SAMPLES;i++)
    {
      for (j = 0;j < N_DIM;j++)
        {
          data->me[j][i] = LENGTH*unirandp();
        }
      fprintf(f, "%f %f %f\n", 
              data->me[0][i], data->me[1][i], data->me[2][i]);
    }
#else

  mu1 = v_get(N_DIM);
  assert(NULL != mu1);
  ret = vec_init(mu1, mu1_dat, N_DIM);
  assert(-1 != ret);
  sigma1 = v_get(N_DIM);
  assert(NULL != sigma1);
  ret = vec_init(sigma1, sigma1_dat, N_DIM);
  assert(-1 != ret);
  rho1 = m_get(N_DIM, N_DIM);
  assert(NULL != rho1);
  ret = mat_init(rho1, rho1_dat, N_DIM*N_DIM);
  assert(-1 != ret);
  ret = rho_check(rho1);
  assert(-1 != ret);

  mu2 = v_get(N_DIM);
  assert(NULL != mu2);
  ret = vec_init(mu2, mu2_dat, N_DIM);
  assert(-1 != ret);
  sigma2 = v_get(N_DIM);
  assert(NULL != sigma2);
  ret = vec_init(sigma2, sigma2_dat, N_DIM);
  assert(-1 != ret);
  rho2 = m_get(N_DIM, N_DIM);
  assert(NULL != rho2);
  ret = mat_init(rho2, rho2_dat, N_DIM*N_DIM);
  assert(-1 != ret);
  ret = rho_check(rho2);
  assert(-1 != ret);

  mu3 = v_get(N_DIM);
  assert(NULL != mu3);
  ret = vec_init(mu3, mu3_dat, N_DIM);
  assert(-1 != ret);
  sigma3 = v_get(N_DIM);
  assert(NULL != sigma3);
  ret = vec_init(sigma3, sigma3_dat, N_DIM);
  assert(-1 != ret);
  rho3 = m_get(N_DIM, N_DIM);
  assert(NULL != rho3);
  ret = mat_init(rho3, rho3_dat, N_DIM*N_DIM);
  assert(-1 != ret);
  ret = rho_check(rho3);
  assert(-1 != ret);

  x = v_get(N_DIM);
  assert(NULL != x);

  for (i = 0;i < N_SAMPLES/3;i++)
    {
      ret = mvnsamp(N_DIM, mu1, sigma1, rho1, x);
      assert(-1 != ret);

      data->me[0][3*i] = x->ve[0];
      data->me[1][3*i] = x->ve[1];
      data->me[2][3*i] = x->ve[2];

      ret = mvnsamp(N_DIM, mu2, sigma2, rho2, x);
      assert(-1 != ret);

      data->me[0][3*i+1] = x->ve[0];
      data->me[1][3*i+1] = x->ve[1];
      data->me[2][3*i+1] = x->ve[2];

      ret = mvnsamp(N_DIM, mu3, sigma3, rho3, x);
      assert(-1 != ret);

      data->me[0][3*i+2] = x->ve[0];
      data->me[1][3*i+2] = x->ve[1];
      data->me[2][3*i+2] = x->ve[2];

      fprintf(f, "%f %f %f\n", 
              data->me[0][3*i], data->me[1][3*i], data->me[2][3*i]);
      fprintf(f, "%f %f %f\n",
              data->me[0][3*i+1], data->me[1][3*i+1], data->me[2][3*i+1]);
      fprintf(f, "%f %f %f\n",
              data->me[0][3*i+2], data->me[1][3*i+2], data->me[2][3*i+2]);
    }

  V_FREE(x);

  M_FREE(rho1);
  V_FREE(sigma1);
  V_FREE(mu1);

  M_FREE(rho2);
  V_FREE(sigma2);
  V_FREE(mu2);

  M_FREE(rho3);
  V_FREE(sigma3);
  V_FREE(mu3);
#endif

  (void)fclose(f);

  mix = mvmixcreate(N_DIM,
                    ranges,
                    N_SAMPLES,
                    data,
                    N_CLASSES, 
                    classes, 
                    init_weights);

  mvmixrand(mix);

  mvmixredraw(mix, DATFILE);
  getchar();

  for (i = 0;i < 100;i++)
    {
      mvmixiter(mix);
      
      mvmixredraw(mix, DATFILE);
      mvmixprint(mix);
      getchar();
    }
  
  mvmixdestroy(mix);
#undef N_DIM
#undef UNIFORM_DISTR
#undef N_SAMPLES
#undef N_CLASSES
}

#undef DATFILE

int main(int argc, char **argv)
{
  int test_nb;

  if (argc != 2)
    {
      fprintf(stderr, "usage: prob_test nb\n");
      exit(1);
    }

  test_nb = atoi(argv[1]);
  
  gopen();

  if (20 == test_nb) test_mvmix();
  if (21 == test_nb) test_mvmix2();

  //gclose();

  //while (1)
  //pause();

  return 0;
}
