/* gfx.c */

#include "mm.h"

static FILE *gpipe = NULL;
static FILE *gpipe2 = NULL;

void gopen()
{
  gpipe = popen("gnuplot", "w");
  //gpipe = stdout;
  assert(NULL != gpipe);
}

int gdraw(char *fmt, ...)
{
  va_list       args;
  int           status;

  va_start(args, fmt);
  status = vfprintf(gpipe, fmt, args);
  va_end(args);

  return status;
}

void gflush()
{
  fflush(gpipe);
}

void gclose()
{
  fclose(gpipe);
  gpipe = NULL;
}

void gopen2()
{
  gpipe2 = popen("gnuplot", "w");
  assert(NULL != gpipe2);
}

int gdraw2(char *fmt, ...)
{
  va_list       args;
  int           status;

  va_start(args, fmt);
  status = vfprintf(gpipe2, fmt, args);
  va_end(args);

  return status;
}

void gflush2()
{
  fflush(gpipe2);
}

void gclose2()
{
  fclose(gpipe2);
  gpipe2 = NULL;
}

