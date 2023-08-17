#include "trig.h"

/**
 * rxengine
 * Calculate the temporal distance from reference point to each point in the field
 * N: number of points in field
 * c: speed of sound [m/s]
 * ref: (x, y, z) coordinate of reference point [m] (3)
 * points: pointer to a matrix of points (N by 3)
 */
float * rxengine(int N, float c, float ** pref, float *** ppoints) {
    // define the output array of tau
    float xdiff, ydiff, zdiff;
    float * tau = (float *) malloc(sizeof(float) * N);
    float * ref = *pref;
    float ** points = *ppoints;

    // iterate through each point
    for(int i = 0; i < N; ++i) {
        xdiff = points[i][0] - ref[0];
        ydiff = points[i][1] - ref[1];
        zdiff = points[i][2] - ref[2];

        tau[i] = sqrtf(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)/c;

        if (PYUSEL_TRIG_DEBUG) {
            printf("%05d: norm(%0.03e, %0.03e, %0.03e)/c = %0.03e us\n", i, xdiff, ydiff, zdiff, 1e6f*tau[i]);
        }
    }

    return tau;
}

/**
 * pwtxengine
 * Calculate the temporal distance from reference point to each point in the field
 * N: number of points in field
 * c: speed of sound [m/s]
 * tref: delay tab of reference element
 * ref: (x, y, z) coordinate of reference point [m]
 * norm: (x, y, z) normal vector
 * xfield, nx, yfield, ny, zfield, nz: pointers to length N arrays of (x, y, z) coordinates in field
 */
float * pwtxengine(int N, float c, float tref, float ** pref, float ** pnorm, float *** ppoints) {
    // iterate through each point
    float xdiff, ydiff, zdiff;
    float * tau = (float *) malloc(sizeof(float) * N);
    float * ref = *pref;
    float * norm = *pnorm;
    float ** points = *ppoints;

    for(int i = 0; i < N; ++i) {
        xdiff = norm[0] * (points[i][0] - ref[0]);
        ydiff = norm[1] * (points[i][1] - ref[1]);
        zdiff = norm[2] * (points[i][2] - ref[2]);

        tau[i] = (xdiff + ydiff + zdiff)/c + tref;

        if (PYUSEL_TRIG_DEBUG) {
            printf("%05d: dot(%0.03e, %0.03e, %0.03e)/c = %0.03e us\n", i, xdiff, ydiff, zdiff, 1e6f*tau[i]);
        }
    }

    return tau;
}

/**
 * genmask3D_new: generate a binary mask of ones or zeros
 * 
 * Parameters:
 * N: number of points
 * fnum: the fnumber along the aperture axis defined by the vector in pnap
 * dyn: whether to use dynamic focussing along the aperture axis
 * pnap: pointer to vector representing the axis of the aperture defined by fnum and dyn
 * pfocus: pointer to vector of length 3 representing the focal point
 * pref: pointer to vector of length 3 representing the reference point
 * ppoints: pointer to N by 3 matrix storing all points
 */
int * genmask3D(int N, float fnum, int dyn, float ** pnap, float ** pfocus, float ** pref, float *** ppoints) {
    float r;
    int in;

    // dereference all vectors
    float * focus = *pfocus;
    float * ref = *pref;
    float * nap = *pnap;
    float ** points = *ppoints;

    int * mask = (int *) malloc(sizeof(int) * N);

    for (int i = 0; i < N; ++i) {
        // calculate radius from center line
        r = nap[0] * (points[i][0] - ref[0]) + nap[1] * (points[i][1] - ref[1]) + nap[2] * (points[i][2] - ref[2]);
        if (r < 0.0f) {r = -r;}

        // determine if within major axis
        in = 0;
        if(0 != dyn) {
            if (2.0f*r <= (points[i][2] - ref[2])/fnum) {in=1;}
        } else {
            if (2.0f*r <= (focus[2] - ref[2])/fnum) {in=1;}
        }

        mask[i] = in;
    }

    return mask;
}
