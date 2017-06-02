
/**
 * @file   prob_test.h
 * @author vr <vr@imagesystems.com>
 * @date   Tue Nov  4 13:14:23 2008
 *
 * @brief  test functions
 *
 *
 */

void gopen(void);
int gdraw(char *fmt, ...) __attribute__((format(printf, 1, 2)));
void gflush(void);
void gclose(void);
void prob_test(void);
void gopen2(void);
int gdraw2(char *fmt, ...) __attribute__((format(printf, 1, 2)));
void gflush2(void);
void gclose2(void);
void prob_test(void);
