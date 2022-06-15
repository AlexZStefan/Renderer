#pragma once
#include "geometry.h"

struct Ray {
    Ray() {};
    Ray(const Point3f& o, const Vec3f& d) : o(o), d(d) {}

    Point3f origin() const { return o; };
    Vec3f direction() const { return d; }

    Point3f at(double t) const {
        return o + t * d;
    }

    Point3f o;
    Vec3f d;
};