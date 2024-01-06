#ifndef __MACHDEP_H__
#define __MACHDEP_H__

//#include <values.h>
#include <float.h>
#include <limits.h>
#define MAXINT INT_MAX

#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

typedef double Real;
#ifdef OLD_COMPLEX
#define Complex complex
#else
#define Complex complex<double>
#endif
#define PATHSEP '/'
#include <stdlib.h>
#include <unistd.h>

#include "fixagg.h"

using namespace std;

#endif
