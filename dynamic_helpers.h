/**
 * This is an implementation of the logistics function, as used
 * by the DF DOP generator.
 *
 * Author: Adrian Tudor Panescu, adrian.panescu@epfl.ch
 */

#ifndef DYNAMIC_HELPERS_H_
#define DYNAMIC_HELPERS_H_

/**
 * Calculates the step that needs to be applied at the current dynamic change
 */
double waveStep(double x, double min, double max, double Y, double scale) {
	return ((x - min) / (max - min)) + Y * scale;
}

/**
 * Calculates the value of the logistics function
 */
double logisticsFunction(double A, double oldY) {
	return A * oldY * (1 - oldY);
}

#endif /* DYNAMIC_HELPERS_H_ */
