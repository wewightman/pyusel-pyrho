#include "rho.h"

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
    for(m=0; m<(M-lag); ++m)
    {
        // get the signals of the elements m and m+lag
        vec1 = (*input)[m];

        // initialize counting variable
        cross = 0;

        // Calculate the mean
        for (n=0; n<N; ++n) cross += vec1[n];
        cross /= (float) N;

        // subtract the mean
        for (n=0; n<N; ++n) vec1[n] -= cross;
    }

    // Calculate and sum the normalized cross correlation for all element pairs
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
    }

    // convert sum to mean
    output /= (float)(M-lag);
    return output;
}