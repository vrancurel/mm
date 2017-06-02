/* mvutils.h */

void mat_print(MAT *mat);
void perm_print(PERM *perm);
void vec_print(VEC *vec);
void vec_print_scilab(VEC *vec, FILE *f);
void mat_print_scilab(MAT *mat, FILE *f);
void vec_parse_scilab(VEC *vec, FILE *f);
void mat_parse_scilab(MAT *mat, FILE *f);
int mat_init(MAT *mat, double *dat, u_int n);
int vec_init(VEC *vec, double *dat, u_int n);
int rho_check(MAT *rho);
void comp_cov(u_int n_dim, VEC *sigma, MAT *rho, MAT *cov);
void comp_rho(u_int n_dim, MAT *cov, VEC *sigma, MAT *rho);
void chol_decomp(MAT *cov);
int nr_fact(int n);
void gen_perm(int k, PERM *sigma);
double diag_prod(MAT *mat);
int comp_det(u_int n_dim, MAT *mat, double *result);
