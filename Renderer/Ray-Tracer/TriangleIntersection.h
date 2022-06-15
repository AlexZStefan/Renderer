#pragma once

#include "Hittable.h"

class TriangleIntersection:public Hittable {
public:
	TriangleIntersection() {};
	//~TriangleIntersection() {};
	TriangleIntersection(Point3f vert0, Point3f vert1, Point3f vert2, shared_ptr<material> m) : v0(vert0), v1(vert1), v2(vert2), mat_ptr(m) {}
	TriangleIntersection(Point3f _v0, Point3f _v1, Point3f _v2, Vec3f _vn0, Vec3f _vn1, Vec3f _vn2, shared_ptr<material> _mat)
		: v0(_v0), v1(_v1), v2(_v2), mat_ptr(_mat), vn0(_vn0), vn1(_vn1), vn2(_vn2) {
		//normal = _vn0 + _vn1 + _vn2 / 3;
	}

	virtual bool bounding_box(AABB& output_box) const override;

	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

	


public:
	Point3f v0;
	Point3f v1;
	Point3f v2;
	Vec3f vn0;
	Vec3f vn1;
	Vec3f vn2;
	Vec3f normal;
	

	shared_ptr<material> mat_ptr;
};


	// compute bounding box points for triangle
inline bool TriangleIntersection::bounding_box(AABB& output_box) const {
	float min[3];
	float max[3];
	for (int i = 0; i < 3; i++) {
		min[i] = std::min(v0[i], std::min(v1[i], v2[i]));
		max[i] = std::max(v0[i], std::max(v1[i], v2[i]));
	}
	output_box = AABB(Vec3f(min[0], min[1], min[2]), Vec3f(max[0], max[1], max[2]));
	return true;
}

// Moller Trumbore Ray Triangle Intersection
// using Crammer`s rule - three vectors (where t, u, v replaces colums to find determinant of vec.x, vec.y, vec.z) / determinant =
// vector C where C = t u v 
// for front faced triangles only - more efficient( / is done at the end since most rays will miss this will be faster)
// only if the determinant is posibive
bool TriangleIntersection::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
	float kEpsilon = 0.00001;
	
	// edges of the triangle
	Vec3f edge1 = v1 - v0;
	Vec3f edge2 = v2 - v0;
	
	// prepare vector to calculate determinant - associated with edge 2
	Vec3f pvec = r.direction().crossProduct(edge2);


	// calculate determinant
	float det = pvec.dotProduct(edge1); // adj

	// if determinant is negative the triangle is backfacing or parallel to ray
	if (det < kEpsilon) return false;
	float invDet = 1 / det;

	Vec3f tvec = r.origin() - v0;

	// u, v = find the Crammer determinant with the scalar triple product by deviding the u determinant by base determinant
	float u = tvec.dotProduct(pvec)* invDet;
	// ray does not intersect triangle - compared to det to avoid the det division - speed up algorithm
	// if u is higher then the det- means that u is higher then 1 
	if (u < 0 || u > 1) return false; 

	Vec3f qvec = tvec.crossProduct(edge1);

	float v = r.direction().dotProduct(qvec) * invDet;
	if (v < 0 || u + v > 1) return false;

	// t is distance from ray origin to the intersection P 
	float t = edge2.dotProduct(qvec) * invDet;

	if (t < 0) return false;

	rec.p = r.at(t);
	rec.t = t;
	
	//rec.normal = (1.0f - u - v) * vn0 + vn1 * u + vn2 * v;
	rec.normal = (vn1 * u) + (vn2 * v) + (vn0 * (1.0f - u - v));
	rec.mat_ptr = mat_ptr;
	
	return true;

}


//bool TriangleIntersection::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
//	// edges of the triangle
//	Vec3f edge2 = v2 - v0;
//	Vec3f edge1 = v1 - v0;
//	
//	// prepare vector to calculate determinant - associated with edge 2
//	Vec3f pvec = r.direction().crossProduct(edge2);
//
//	// calculate determinant
//	float det = edge1.dotProduct(pvec); // adj
//
//	float kEpsilon = 0.00001;
//
//	// if determinant is negative the triangle is backfacing or parallel to ray
//	if (det < kEpsilon) return false; 
//
//	float invDet = 1 / det; 
//
//	Vec3f tvec = r.origin() - v0;
//	// u, v = find the Crammer determinant with the scalar triple product by deviding the u determinant by base determinant
//	float u = tvec.dotProduct(pvec) * invDet; 
//	if (u < 0 || u > 1) return false; // ray does not intersect triangle
//	
//	Vec3f qvec = tvec.crossProduct(edge1);
//
//	v = r.direction().dotProduct(qvec) * invDet;
//
//	// t is distance from ray origin to the intersection P 
//	float t = edge2.dotProduct(pvec) * invDet;
//
//}
