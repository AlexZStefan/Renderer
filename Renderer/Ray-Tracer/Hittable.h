#pragma once
#include "Ray.h"
#include "common.h"
#include "AABB.h"

class material;

struct hit_record {
	Point3f p;
	Vec3f normal;
	double t;
	bool front_face;
	double u, v; 

	shared_ptr<material> mat_ptr;

	inline void set_face_normal(const Ray& r, const Vec3f& outward_normal) {
		// depending on hit face -inner our outer - returns normal
		front_face = (r.direction().dotProduct(outward_normal)) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

class Hittable
{
public:
	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const = 0;
	virtual bool bounding_box(AABB& output_box) const = 0;
};

