/* mvmoduni.c */

#include "mm.h"
#include <math.h>

static int densityfunc(tmvmix *mix,
                       tmvmodinst *model, 
                       void *unused,
                       VEC *sample,
                       double *result)
{
  double prod;
  u_int i;

  prod = 1;
  for (i = 0;i < mix->n_dim;i++)
    prod *= mix->ranges[i].end - mix->ranges[i].start;

  *result = 1.0/prod;

  return 0;
}


static int jointcdffunc(tmvmix *mix,
                        tmvmodinst *model, 
                        void *ddp,
                        VEC *sample, 
                        double *result)
{
  double prod;
  u_int i;

  prod = 1;
  for (i = 0; i < mix->n_dim; ++i)
    prod  *= (sample->ve[i] - mix->ranges[i].start) /
      (mix->ranges[i].end - mix->ranges[i].start);

  *result = prod;
  return 0;
}

tmvmodclass mvmoduniclass = 
  {
    .name = "uniform",
    .size = 0,
    .init = NULL,
    .reset = NULL,
    .copy = NULL,
    .rand = NULL,
    .densityinit = NULL,
    .density = densityfunc,
    .jointcdf = jointcdffunc,
    .densityfree = NULL,
    .update = NULL,
    .print = NULL,
    .redraw = NULL,
    .save = NULL,
    .restore = NULL,
    .destroy = NULL,
  };
