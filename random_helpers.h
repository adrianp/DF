/**
 * Some useful functions for generating random stuff.
 *
 * Author: Adrian Tudor Panescu, adrian.panescu@epfl.ch
 */

#include <stdlib.h>

#ifndef RANDOM_HELPERS_H_
#define RANDOM_HELPERS_H_

/**
 * Generates a random double
 */
double random_double() {
	return rand() / double(RAND_MAX);
}

/**
 * Generates a random double in the interval [min, max]
 */
double random_double(double min, double max) {
	return (max - min) * random_double() + min;
}

/**
 * Generates a random double in the interval [min, max)
 */
double random_double_exclusive(double min, double max) {
	double res = (max - min) * random_double() + min;
	while (fabs(max - res) < std::numeric_limits<double>::epsilon()) {
		res = (max - min) * random_double() + min;
	}
	return res;
}

// TODO (adrianp): http://www.cplusplus.com/reference/algorithm/random_shuffle/

#endif // RANDOM_HELPERS_H_
