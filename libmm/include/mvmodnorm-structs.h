/* mvmodnorm-structs.h */

typedef struct
{
  VEC                   *mu;			/*!< means [n_dim] */
  MAT                   *cov;         /*!< covariance matrix [n_dim][n_dim] */
} tmvmodnormparams;
