/* mvmix-structs.h */

typedef struct
{
  double        start;
  double        end;
} tmvrange;

struct smvmix;

typedef struct
{
  double			weight;		/*!< weight */
  double			*probs;	      /*!< probabilities [n_samples] */
} tmvmodinst;

//int functions returns 0 if OK, else -1

typedef int (*tmvmodinitfunc)(struct smvmix *mix, tmvmodinst *model);

typedef void (*tmvmodresetfunc)(struct smvmix *mix, tmvmodinst *model);

typedef void (*tmvmodcopyfunc)(struct smvmix *mix, tmvmodinst *new, tmvmodinst *old);

typedef void (*tmvmodrandfunc)(struct smvmix *mix, tmvmodinst *model);

typedef int (*tmvmoddensityinitfunc)(struct smvmix *mix,
                                     tmvmodinst *model, 
                                     void **ddp); //density data

typedef int (*tmvmoddensityfunc)(struct smvmix *mix,
                                 tmvmodinst *model, 
                                 void *dd, //density data
                                 VEC *sample,
                                 double *result);

typedef int (*tmvmodjointcdffunc)(struct smvmix *mvmix,
                                  tmvmodinst *model, 
                                  void *dd, //density data
                                  VEC *sample,
                                  double *result);

typedef void (*tmvmoddensityfreefunc)(struct smvmix *mix,
                                      tmvmodinst *model, 
                                      void *dd); //density data

typedef int (*tmvmodupdatefunc)(struct smvmix *mix,
                                tmvmodinst *model);

typedef void (*tmvmodprintfunc)(struct smvmix *mix,
                                tmvmodinst *model);

typedef void (*tmvmodredrawfunc)(struct smvmix *mix,
                                 tmvmodinst *model);

typedef void (*tmvmodsavefunc)(struct smvmix *mix,
                               tmvmodinst *model,
                               FILE *f);

typedef void (*tmvmodrestorefunc)(struct smvmix *mix,
                                  tmvmodinst *model,
                                  FILE *f);

typedef void (*tmvmoddestroyfunc)(struct smvmix *mix,
                                  tmvmodinst *model);

typedef struct
{
  char                          *name;
  u_int                         size;           /*!< size of instance */
  tmvmodinitfunc		init;		/*!< create an instance */
  tmvmodresetfunc		reset;		/*!< reset parameters */
  tmvmodcopyfunc		copy;		/*!< copy parameters */
  tmvmodrandfunc		rand;           /*!< set random parameters */
  tmvmoddensityinitfunc		densityinit;	/*!< init density */
  tmvmoddensityfunc		density;	/*!< compute density */
  tmvmodjointcdffunc            jointcdf;       /*!< compute joint cdf */
  tmvmoddensityfreefunc		densityfree;	/*!< free density */
  tmvmodupdatefunc		update;    /*!< update params with new probs */
  tmvmodprintfunc               print;     /*!< print debug information */
  tmvmodredrawfunc              redraw;     /*!< redraw parameters (gnuplot) */
  tmvmodsavefunc                save;           /*!< save parameters */
  tmvmodrestorefunc             restore;        /*!< restore parameters */
  tmvmoddestroyfunc		destroy;	/*!< destroy instance */
} tmvmodclass;

typedef struct smvmix
{
  u_int                         n_dim;          /*!< number of dimensions */
  const tmvrange                *ranges;  /*!< hyperspace definition [n_dim] */
  u_int				n_samples;	/*!< number of samples */
  const MAT			*data;          /*!< data [n_dim][n_samples] */
  u_int				n_kernels;	/*!< number of kernels */
  const tmvmodclass		**classes;	/*!< classes [n_kernels] */
  const double		     *rel_weights; /*!< relative weights [n_kernels] */
  tmvmodinst			**insts;	/*!< instances [n_kernels] */
  double			*probs;	      /*!< probabilities [n_samples] */
} tmvmix;

typedef struct
{
  tmvmodinst                    **insts;
  double                        *probs;
} tmvmixsnapshot;
