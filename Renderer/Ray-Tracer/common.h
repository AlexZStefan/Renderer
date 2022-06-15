#pragma once

#define _USE_MATH_DEFINES
#include "Ray.h"
#include "geometry.h"

#include <cstdlib>
#include <memory>
#include <limits>
#include <cmath>
#include <math.h>

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

const double infinity = std::numeric_limits<double>::infinity();
const double pi = M_PI;

inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0f;
}
// (return a random < 1) + 1;
inline double random_double() {
	return rand() / (RAND_MAX + 1.0f);
}

inline double random_double(double min, double max) {
	return min + (max - min) * random_double();
}

