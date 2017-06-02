/* mvmix.h */

tmvmix *mvmixcreate(const u_int n_dim,
                    const tmvrange *ranges,
                    const u_int n_samples,
                    const MAT *data,
                    const u_int n_kernels,
                    const tmvmodclass **classes,
                    const double *init_weights);

void mvmixreset(tmvmix *mix);

void mvmixrand(tmvmix *mix);

int mvmixprob(tmvmix *mix);

int mvmixiter(tmvmix *mix);

void mvmixprint(tmvmix *mix);

void mvmixredraw(tmvmix *mix, const char *datfile);

void mvmixdestroy(tmvmix *mix);

double mvmixlikelihood(tmvmix *mix);

tmvmixsnapshot *mvmixsnapshotcreate(tmvmix *mvmix);

void mvmixsnapshotsave(tmvmix *mix, tmvmixsnapshot *snapshot);

void mvmixsnapshotrestore(tmvmix *mix, tmvmixsnapshot *snapshot);

void mvmixsnapshotdestroy(tmvmix *mix, tmvmixsnapshot *snapshot);

int mvmixsave(tmvmix *mix, char *ofname);

int mvmixrestore(tmvmix *mix, char *ifname);
