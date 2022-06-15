#pragma once

#include "Material.h"

class Sphere : public Hittable {
public:
	Sphere() {};

	//Sphere(const Point3f& c, double r) : centre(c), radius(r) {}

	Sphere(Point3f cen, double r, shared_ptr<material> m) : centre(cen), radius(r), mat_ptr(m){}

	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

	virtual bool bounding_box(AABB& output_box) const override;

	Point3f centre;
	double radius;
	std::shared_ptr<material> mat_ptr;

private:
	// compute anglefor a given point on the unit sphere centred at the origin
	static void get_sphere_uv(const Point3f& p, double& u, double& v) {
		
			// p - point on the sphere of radius one, centered at the origin.
			// u - returned value [0,1] of angle around the Y axis from X=-1.
			// v - returned value [0,1] of angle from Y=-1 to Y=+1.
	
		auto theta = acos(-p.y);
		auto phi = atan2(-p.z, p.x) + M_PI;
		// inherited from hittable 
		u = phi / (2 * M_PI);
		v = theta / M_PI;
	}

};
	
inline bool Sphere::bounding_box(AABB& output_box) const {
	// create the 2 points that make up the surrounding box
	output_box = AABB(centre - Vec3f(radius, radius, radius),
		centre + Vec3f(radius, radius, radius));
	return true;
}


bool Sphere::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
	// ray origin to centre of sphere
	Vec3f oc = r.origin() - centre;
	// 
	auto a = r.direction().norm();
	// is 
	auto half_b = oc.dotProduct(r.direction());
	
	// centre of sphere 
	auto c = oc.norm() - radius * radius;

	auto discriminant = half_b * half_b - a * c; // a=1 as ray is normalised​

	if (discriminant < 0) return false;

	auto sqrtd = sqrt(discriminant);

	// Find nearest root that lies in the acceptable range
	// where ray is at a specific time
	auto root = (-half_b - sqrtd) / a;

	if (root < t_min || t_max < root) {
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root)
			return false;
	};

	// update hit record data 
	rec.t = root;
	rec.p = r.at(rec.t);
	// Calculate normal from where point hits the sphere - the centre of sphere
	Vec3f outward_normal = (rec.p - centre) / radius;
	rec.set_face_normal(r, outward_normal);
	// get sphere uv 
	get_sphere_uv(outward_normal, rec.u, rec.v);
	rec.mat_ptr = mat_ptr;
	
	return true;
}