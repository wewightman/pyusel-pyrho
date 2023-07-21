#ifndef ___pyusel_pyrho___
#define ___pyusel_pyrho___
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern float lagNRho(float *** input, int M, int N, int lag);

#endif

#ifndef PYUSEL_DEBUG
#define PYUSEL_DEBUG 0
#endif