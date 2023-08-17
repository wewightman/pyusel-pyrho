#ifndef ___pyusel_pyrho_trig___
#define ___pyusel_pyrho_trig___
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern float * rxengine(int N, float c, float ** pref, float *** ppoints);
extern float * pwtxengine(int N, float c, float ** pref, float ** pnorm, float *** ppoints);
extern int * genmask3D(int N, float fnum, int dyn, float ** pnap, float ** pfocus, float ** pref, float *** ppoints);

#endif

#ifndef PYUSEL_TRIG_DEBUG
#define PYUSEL_TRIG_DEBUG 0
#endif