#include "rho.h"
#include "cubic.h"

/**
 * lagnrho: calculate the nth lag
 * 
 * Parameters:
 * input: pointer to pointer-of-pointers-of-floats (2D matrix) of size M by N. m = channel index, n = sample index
 * M: the number of vectors in the matrix
 * N: the number of samples in each vector
 * lag: the lag index to form pairs over
 * 
 * Returns:
 *  output: a float containing the lag(lag) coherence
*/
float lagNRho(float *** input, int M, int N, int lag)
{
    float * vec1;
    float * vec2;
    float cross, self1, self2, output;
    int m, n;

    // calculate and subtract the mean of each element
    if (PYUSEL_DEBUG) printf("Subtracting avg:\n");
    for(m=0; m<(M-lag); ++m)
    {
        // get the signals of the elements m and m+lag
        vec1 = (*input)[m];

        // initialize counting variable
        cross = 0;

        // Calculate the mean
        for (n=0; n<N; ++n) cross += vec1[n];
        if (PYUSEL_DEBUG) printf("  vector sum: %e\n", cross);
        cross /= (float) N;
        if (PYUSEL_DEBUG) printf("  vector avg: %e\n", cross);

        // subtract the mean
        for (n=0; n<N; ++n) vec1[n] -= cross;
    }

    // Calculate and sum the normalized cross correlation for all element pairs
    if (PYUSEL_DEBUG) printf("Calculating norm x-corr:\n");
    output = 0.0f;
    for(m=0; m<(M-lag); ++m)
    {
        // get the signals of the elements m and m+lag
        vec1 = (*input)[m];
        vec2 = (*input)[m+lag];
        
        // clear the cross and self correlation terms for this element pair
        cross = 0.0f;
        self1 = 0.0f;
        self2 = 0.0f;

        // for each sample
        for (n=0; n<N; ++n)
        {
            self1 += vec1[n] * vec1[n]; // variance of vec1
            cross += vec1[n] * vec2[n]; // covariance of vec1 x vec2
            self2 += vec2[n] * vec2[n]; // variance of vec2
        }

        // calculate normalized cross correllation and sum with previous
        output += cross/sqrtf(self1*self2);
        if (PYUSEL_DEBUG) printf("  This covariance is: %e\n", cross/sqrtf(self1*self2));
    }

    if (PYUSEL_DEBUG) printf("Pre-division: %e\n", output);

    // convert sum to mean
    output /= (float)(M-lag);

    if (PYUSEL_DEBUG) printf("Post division: %e\n", output);
    if (PYUSEL_DEBUG) printf("Post division: %p\n", (void*) &output);
    return output;
}

void summat(int M, int N, float ** mat1, float ** mat2, float ** out);
float ** getmat(int M, int N);
void freemat(int M, float ** mat);
float *** gettsel(int Ntx, int Nrx, int ip, int Ntsel, float * tsel, float ** tautx, float ** taurx);
void freetsel(int Ntx, int Nrx, float *** tsel);

/**
 * lagNRhoSet: convert given transmit and recieve delay tabs for each transmit, each receive, and each point, compute lagN coherence
 * 
 * Parameters:
 * lag: lag index to use
 * Np: number of points
 * Ntx: number of transmissions
 * Nrx: number of receives
 * Nt: number of time points
 * Ts: sampling period
 * tstart: starting value for time traces
 * tautx: pointer to Ntx by Np matrix of transmist delays
 * taurx: pointer to Nrx by Np matrix of receive delays
 * rf: data tensor  Ntx by Nrx by Nt
 * 
 * Returns:
 * the lagN vlaue for each interpolator
 */
float * lagNRhoSet(int lag, int Np, int Ntx, int Nrx, int Nt, float Ts, float tstart, float dtw, float *** tautx, float *** taurx, float **** rf)
{

    if (PYUSEL_DEBUG || PYUSEL_PYRHO_DEBUG) printf("Started initialization...");
    int ip, itx, irx, it;
    float * lags = (float*) malloc(sizeof(float) * Np);


    int Ntsel = ((int) (dtw/Ts + 1));
    float * tsel = (float *) malloc(sizeof(float) * (1+2*Ntsel));

    float ** extracted = (float **) malloc(sizeof(float*) * Nrx);
    
    float *** tsels;
    float ** buffer;
    
    dtw = Ts * (float)Ntsel;

    if (PYUSEL_DEBUG || PYUSEL_PYRHO_DEBUG) printf("Finished initializing");

    // make interpolators
    IntrpData1D_Fixed *** interps = (IntrpData1D_Fixed ***) malloc(sizeof(IntrpData1D_Fixed **) * Ntx);
    for (itx = 0; itx < Ntx; ++itx)
    {
        interps[itx] = (IntrpData1D_Fixed **) malloc(sizeof(IntrpData1D_Fixed *) * Nrx);
        for (irx = 0; irx < Nrx; ++irx)
        {
            interps[itx][irx] = tie_knots1D_fixed(&((*rf)[itx][irx]), Nt, Ts, tstart, 0.0f, 1);
        }
    }

    if (PYUSEL_DEBUG || PYUSEL_PYRHO_DEBUG) printf("Finished initializing interpolators");

    // define tsel
    for (it = -Ntsel; it <= Ntsel; ++it) tsel[it] = dtw * (float)it;

    // iterate through each point
    if (PYUSEL_DEBUG || PYUSEL_PYRHO_DEBUG) printf("Iterating through points");
    for (ip = 0; ip < Np; ++ip)
    {
        tsels = gettsel(Ntx, Nrx, ip, Ntsel*2+1, tsel, *tautx, *taurx);
        buffer = (float **) malloc(sizeof(float*) * Nrx);
        for (itx = 0; itx < Ntx; ++itx) buffer[itx] = (float*) calloc(1+2*Ntsel, sizeof(float));
        
        if (PYUSEL_DEBUG || PYUSEL_PYRHO_DEBUG) printf("  %08d: summing transmissions together", ip);
        // sum all transmission events together
        for (itx = 0; itx < Ntx; ++itx)
        {
            for (irx = 0; irx < Nrx; ++irx) extracted[irx] = cubic1D_fixed(&(tsels[itx][irx]), Ntsel*2+1, &(interps[itx][irx]));

            summat(irx, 1+2*Ntsel, buffer, extracted, buffer);

            for (irx = 0; irx < Nrx; ++irx) free(extracted[irx]);
        }

        if (PYUSEL_DEBUG || PYUSEL_PYRHO_DEBUG) printf("    Calculating LNC for this point");
        lags[ip] = lagNRho(&(buffer), Nrx, 1+2*Ntsel, lag);

        freetsel(Ntx, Nrx, tsels);
        freemat(Nrx, buffer);
    }

    // clean up all intermediate buffers
    free(extracted);
    free(tsel);

    return lags;
}

void summat(int M, int N, float ** mat1, float ** mat2, float ** out)
{
    int m, n;
    // sum matrix
    for (m = 0; m < M; ++m) 
    {
        for (n=0; n<N; ++n) 
        {
            out[m][n] = mat1[m][n] + mat2[m][n];
        }
    }
}

float ** getmat(int M, int N)
{
    int m;
    // initialize an empty matrix
    float ** out = (float **) malloc(sizeof(float *) * M);
    for(m = 0; m < M; ++m) out[m] = (float *) malloc(sizeof(float) * N);
}

void freemat(int M, float ** mat)
{
    // initialize an empty matrix
    for(int m = 0; m < M; ++m) free(mat[m]);
    free(mat);
}

float *** gettsel(int Ntx, int Nrx, int ip, int Ntsel, float * tsel, float ** tautx, float ** taurx)
{
    int itx, irx, it;
    float tx0, rx0;
    float *** tsels  = (float ***) malloc(sizeof(float **) * Ntx);
    for (itx = 0; itx < Ntx; itx++) tsels[itx] = getmat(Nrx, Ntsel);

    for (itx = 0; itx < Ntx; ++itx)
    {
        tx0 = tautx[itx][ip];
        for (irx = 0; irx < Nrx; ++irx)
        {
            rx0 = taurx[itx][ip];
            for (it = 0; it < Ntsel; ++it)
            {
                tsels[itx][irx][it] = tx0 + rx0 + tsel[it];
            }
        }
    }

    return tsels;
}

void freetsel(int Ntx, int Nrx, float *** tsel)
{
    int itx, irx;

    for (itx = 0; itx < Ntx; ++itx)
    {
        for (irx = 0; irx < Nrx; ++irx)
        {
            free(tsel[itx][irx]);
        }
        free(tsel[itx]);
    }
    free(tsel);
}