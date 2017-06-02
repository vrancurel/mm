/* mvmix.c */

#include "mm.h"
#include <math.h>

/** 
 * create a mixture
 * 
 * @param n_kernels 
 * @param classes 
 * 
 * @return a newly allocated mixture
 */
tmvmix *mvmixcreate(const u_int n_dim,
                    const tmvrange *ranges,
                    const u_int n_samples,
                    const MAT *data,
                    const u_int n_kernels,
                    const tmvmodclass **classes,
                    const double *rel_weights)
{
  tmvmix *mix;
  u_int i;
  double sum;
  int ret;

  mix = xmalloc(sizeof (*mix));

  memset(mix, 0, sizeof (*mix));

  mix->n_dim = n_dim;
  mix->ranges = ranges;
  mix->n_samples = n_samples;
  mix->data = data;
  mix->n_kernels = n_kernels;
  mix->classes = classes;
  mix->rel_weights = rel_weights;

  mix->insts = xmalloc(n_kernels * sizeof (void *));

  memset(mix->insts, 0, n_kernels * sizeof (void *));

  mix->probs = xmalloc(n_samples * sizeof (double));

  memset(mix->probs, 0, n_samples * sizeof (double));

  sum = 0.0;
  for (i = 0;i < mix->n_kernels;i++)
    {
      mix->insts[i] = xmalloc(sizeof (tmvmodinst) + 
                              mix->classes[i]->size);

      mix->insts[i]->probs = xmalloc(n_samples * sizeof (double));

      memset(mix->insts[i]->probs, 0, n_samples * sizeof (double));

      if (NULL != mix->classes[i]->init)
        {
          ret = mix->classes[i]->init(mix, mix->insts[i]);
          assert(-1 != ret);
        }

      sum += rel_weights[i];
    }

 for (i = 0;i < mix->n_kernels;i++)
   mix->insts[i]->weight = rel_weights[i] / sum;

  return mix;
}

/** 
 * reset all parameters
 * 
 * @param mix 
 */
void mvmixreset(tmvmix *mix)
{
  u_int i;
  double sum;

  memset(mix->probs, 0, mix->n_samples * sizeof (double));

  sum = 0.;
  for (i = 0;i < mix->n_kernels;i++)
    {
      if (NULL != mix->classes[i]->reset)
        mix->classes[i]->reset(mix, mix->insts[i]);

      memset(mix->insts[i]->probs, 0, mix->n_samples * sizeof (double));

      sum += mix->rel_weights[i];
    }

 for (i = 0;i < mix->n_kernels;i++)
   mix->insts[i]->weight = mix->rel_weights[i] / sum;
}

/** 
 * set random parameters
 * 
 * @param mix 
 */
void mvmixrand(tmvmix *mix)
{
  u_int i;

  for (i = 0;i < mix->n_kernels;i++)
    {
      if (NULL != mix->classes[i]->rand)
        mix->classes[i]->rand(mix, mix->insts[i]);
    }
}

int mvmixprob(tmvmix *mix)
{
  u_int i, j;
  int ret;
  
  memset(mix->probs, 0, mix->n_samples * sizeof (double));
  
  for (i = 0;i < mix->n_kernels;i++)
    {
      void *dd;

      memset(mix->insts[i]->probs, 0, mix->n_samples * sizeof (double));

      if (NULL != mix->classes[i]->densityinit)
        {
          ret = mix->classes[i]->densityinit(mix, mix->insts[i], &dd);
          if (-1 == ret)
            {
              fprintf(stderr, "density init failed\n");
              return -1;
            }
        }

      for (j = 0;j < mix->n_samples;j++)
	{
          VEC *sample;
          int k;
          double result;

	  /*
	   * compute probabilities for the data set given the kernel and
	   * its weight
	   */
          sample = v_get(mix->n_dim);
          assert(NULL != sample);

          for (k = 0;k < mix->n_dim;k++)
            sample->ve[k] = mix->data->me[k][j];
          
          ret = mix->classes[i]->density(mix,
                                         mix->insts[i], 
                                         dd,
                                         sample,
                                         &result); 
          assert(-1 != ret);

          mix->insts[i]->probs[j] = mix->insts[i]->weight * result;
          
          V_FREE(sample);

	  /*
	   * compute total probabilities for each sample
	   */
	  mix->probs[j] += mix->insts[i]->probs[j];
	}
      
      if (NULL != mix->classes[i]->densityfree)
        {
          mix->classes[i]->densityfree(mix, mix->insts[i], dd);
        }
    }

  return 0;
}

/** 
 * Perform one iteration
 * 
 * @param mix 
 *
 * @return 0 if OK, else error, e.g. density not computable
 */
int mvmixiter(tmvmix *mix)
{
  u_int i, j;
  double sum;
  int ret;

  ret = mvmixprob(mix);
  if (0 != ret)
    return -1;

  /*
   * adjust probabilities, recompute weight
   */
  sum = 0.0;
  for (i = 0;i < mix->n_kernels;i++)
    {
      mix->insts[i]->weight = 0.0;
      for (j = 0;j < mix->n_samples;j++)
	{
	  mix->insts[i]->probs[j] /= mix->probs[j];
	  mix->insts[i]->weight += mix->insts[i]->probs[j];
	}
      
      mix->insts[i]->weight /= mix->n_samples;

      sum += mix->rel_weights[i];
    }

  /*
   * update parameters
   */
  for (i = 0;i < mix->n_kernels;i++)
    {
      if (NULL != mix->classes[i]->update)
        {
          ret = mix->classes[i]->update(mix, mix->insts[i]);
          assert(-1 != ret);
        }

      mix->insts[i]->weight = 
	0.9 * mix->insts[i]->weight + 
	0.1 * mix->rel_weights[i] / sum;
    }

  return 0;
}

/** 
 * print weights and parameters
 * 
 * @param mix 
 */
void mvmixprint(tmvmix *mix)
{
  u_int i;
  
  for (i = 0;i < mix->n_kernels;i++)
    {
      printf("kernel %d weight %f\n", i, mix->insts[i]->weight);

      if (NULL != mix->classes[i]->print)
        {
          mix->classes[i]->print(mix, mix->insts[i]);
        }
      printf("--\n");
    }

  printf("likelihood %f\n", mvmixlikelihood(mix));
  printf("---\n");
}

/** 
 * redraw density ellipses (only works with 2-D and 3-D)
 * 
 * @param mix 
 * @param datfile 
 */
void mvmixredraw(tmvmix *mix, const char *datfile)
{
  u_int i;
  int special = 0;

  if (!(mix->n_dim == 1 ||
        mix->n_dim == 2 ||
        mix->n_dim == 3))
    {
      printf("unable to draw in such dimension\n");
      return ;
    }

  if (mix->n_dim == 1)
    {
      gdraw("plot ");
      special = 1;
    }
  else if (mix->n_dim == 2)
    gdraw("plot \"%s\"", datfile);
  else if (mix->n_dim == 3)
    gdraw("splot \"%s\"", datfile);

  for (i = 0;i < mix->n_kernels;i++)
    {
      if (NULL != mix->classes[i]->redraw)
        {
          if (mix->n_dim == 1 || !special)
            gdraw(",");
          else
            special = 0;

          mix->classes[i]->redraw(mix, mix->insts[i]);
        }
    }
  gdraw("\n");
  gflush();
}

/** 
 * destroy the mixture
 * 
 * @param mix 
 */
void mvmixdestroy(tmvmix *mix)
{
  u_int i;

  for (i = 0;i < mix->n_kernels;i++)
    {
      free(mix->insts[i]->probs);

      if (NULL != mix->classes[i]->destroy)
        {
          mix->classes[i]->destroy(mix, mix->insts[i]);
        }
    }

  free(mix->probs);
  free(mix->insts);
  free(mix);
}

/** 
 * compute the likelihood. need at least one iteration
 * 
 * @param mix 
 * 
 * @return 
 */
double mvmixlikelihood(tmvmix *mix)
{
  double tmp;
  u_int i;

  tmp = 0.0;
  for (i = 0;i < mix->n_samples;i++)
    tmp += log(mix->probs[i]);

  return tmp / mix->n_samples;
}

/** 
 * create a snapshot
 * 
 * @param mvmix 
 * 
 * @return 
 */
tmvmixsnapshot *mvmixsnapshotcreate(tmvmix *mix)
{
  tmvmixsnapshot *snapshot;
  u_int i;
  int ret;
  
  snapshot = xmalloc(sizeof (*snapshot));
  
  snapshot->insts = xmalloc(mix->n_kernels * sizeof (void *));
  
  snapshot->probs = xmalloc(mix->n_samples * sizeof (double));
  
  for (i = 0;i < mix->n_kernels;i++)
    {
      snapshot->insts[i] = xmalloc(sizeof (tmvmodinst) + 
                                   mix->classes[i]->size);

      snapshot->insts[i]->probs = 
	xmalloc(mix->n_samples * sizeof (double));

      if (NULL != mix->classes[i]->init)
        {
          ret = mix->classes[i]->init(mix, snapshot->insts[i]);
          assert(-1 != ret);
        }
    }

  return snapshot;
}

/** 
 * save the current estimates (weights, parameters and probabilities) in
 * snapshot
 * 
 * @param mix 
 * @param state 
 * 
 * @return 
 */
void mvmixsnapshotsave(tmvmix *mix, tmvmixsnapshot *snapshot)
{
  u_int i;

  memcpy(snapshot->probs, mix->probs, 
         mix->n_samples * sizeof (double));
  
  for (i = 0;i < mix->n_kernels;i++)
    {
      snapshot->insts[i]->weight = mix->insts[i]->weight;

      memcpy(snapshot->insts[i]->probs,
             mix->insts[i]->probs,
             mix->n_samples * sizeof (double));

      if (NULL != mix->classes[i]->copy)
        {
          mix->classes[i]->copy(mix,
                                snapshot->insts[i], 
                                mix->insts[i]);
        }
    }
}

/** 
 * restore the saved estimates
 * 
 * @param mix 
 * @param state 
 * 
 * @return 
 */
void mvmixsnapshotrestore(tmvmix *mix, tmvmixsnapshot *snapshot)
{
  u_int i;
 
  memcpy(mix->probs, snapshot->probs, 
         mix->n_samples * sizeof (double));
  
  for (i = 0;i < mix->n_kernels;i++)
    {
      mix->insts[i]->weight = snapshot->insts[i]->weight;

      memcpy(mix->insts[i]->probs,
             snapshot->insts[i]->probs,
             mix->n_samples * sizeof (double));

      if (NULL != mix->classes[i]->copy)
        {
          mix->classes[i]->copy(mix,
                                mix->insts[i], 
                                snapshot->insts[i]);
        }
    }
}

/** 
 * destroy a snapshot
 * 
 * @param mix 
 * @param snapshot 
 */
void mvmixsnapshotdestroy(tmvmix *mix, tmvmixsnapshot *snapshot)
{
  u_int i;

  for (i = 0;i < mix->n_kernels;i++)
    {
      free(snapshot->insts[i]->probs);

      if (NULL != mix->classes[i]->destroy)
        {
          mix->classes[i]->destroy(mix, snapshot->insts[i]);
        }
    }

  free(snapshot->probs);
  free(snapshot->insts);
  free(snapshot);
}

/*
 *
 */

int mvmixsave(tmvmix *mix, char *ofname)
{
  FILE *ofile;
  u_int i;
  
  ofile = fopen(ofname, "w");
  if (NULL == ofile)
    return -1;

  for (i = 0;i < mix->n_kernels;i++)
    {
      fprintf(ofile, "%d\n", i);
      fprintf(ofile, "%f\n", mix->insts[i]->weight);
      if (NULL != mix->classes[i]->save)
        {
          mix->classes[i]->save(mix, mix->insts[i], ofile);
        }
    }

  fclose(ofile);

  return 0;
}

int mvmixrestore(tmvmix *mix, char *ifname)
{
  FILE *ifile;
  u_int i;
  
  ifile = fopen(ifname, "r");
  if (NULL == ifile)
    return -1;

  for (i = 0;i < mix->n_kernels;i++)
    {
      int nb;

      printf("restoring %d\n", i);

      fscanf(ifile, "%d", &nb);
      if (i != nb)
        {
          printf("%d != %d\n", nb, i);
          exit(1);
        }
      fscanf(ifile, "%lf", &mix->insts[i]->weight);
      if (NULL != mix->classes[i]->restore)
        {
          mix->classes[i]->restore(mix, mix->insts[i], ifile);
        }
      
      printf("ok\n");
    }

  //fclose(ifile);

  return 0;
}
