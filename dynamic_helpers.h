#ifndef DYNAMIC_HELPERS_H_
#define DYNAMIC_HELPERS_H_

double waveStep(double x, double min, double max, double Y, double scale)
{
    return ((x-min) / (max-min)) + Y * scale;
}

double logisticsFunction(double A, double oldY)
{
    return A * oldY * (1 - oldY);
}

#endif /* DYNAMIC_HELPERS_H_ */
