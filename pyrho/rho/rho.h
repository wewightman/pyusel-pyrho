#ifndef ___pyusel_pyrho_rho___
#define ___pyusel_pyrho_rho___
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cubic.h"

extern float * lagNRhoSet(int lag, int Np, int Ntx, int Nrx, int Nt, float Ts, float tstart,  float dtw, float *** tautx, float *** taurx, float **** rf);
extern float lagNRho(float *** input, int M, int N, int lag);

#endif

#ifndef PYUSEL_DEBUG
#define PYUSEL_DEBUG 0
#endif

#ifndef PYUSEL_PYRHO_DEBUG
#define PYUSEL_PYRHO_DEBUG 1
#endif