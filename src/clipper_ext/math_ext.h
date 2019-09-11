#include <iostream>
#include <vector>
// Quick-and-dirty approximation to the inverse error function based on
// "A handy approximation for the error function and its inverse" by
// Sergei Winitzki. Found at
// https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
// Maximum error is about 6e-3.
double approx_inverse_erf(double x)
{
    long double tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0 : 1.0;
    lnx = log(1.0L-x*x);
    tt1 = 2/(M_PI*0.147) + 0.5 * lnx;
    tt2 = 1/0.147 * lnx;
    return (sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
}

double approx_inverse_erfc(double x)
{
    return approx_inverse_erf(1.0-x);
}
