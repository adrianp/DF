#include <stdlib.h>

#ifndef RANDOM_HELPER_H_
#define RANDOM_HELPER_H_

double random_double()
{
    return rand() / double(RAND_MAX);
}

double random_double(double min, double max)
{
    return (max - min) * random_double() + min;
}

double random_double_exclusive(double min, double max)
{
    double res = (max - min) * random_double() + min;
    while(res == max) // TODO: comparing floating point with == is unsafe
    {
        res = (max - min) * random_double() + min;
    }
    return res;
}

#endif // RANDOM_HELPER_H_
