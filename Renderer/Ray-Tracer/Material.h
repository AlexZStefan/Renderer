#pragma once

#include "common.h"  
#include "geometry.h"
#include "Hittable.h"
#include "texture.h"

Vec3f reflect(const Vec3f& v, const Vec3f& n) {
	return v - 2 * v.dotProduct(n) * n;
};


// computed from direction of ray_in.normalized, record_hit_normal and ratio of refraction
Vec3f refract(const Vec3f& uv, const Vec3f& n, double etai_over_etat) {
	// should the R inside sqrt be R`|| ? 
	// R`= q/q` * ( R + cos0) + -sqrt(1-|R` _|_|^2 ; R`_|_ = q/q`(R +cos0)n;
	// calculate the angle above surface
	auto cos_theta = fmin(-uv.dotProduct(n), 1.0);
	// calc R`_|_ 
	Vec3f r_out_perp = etai_over_etat * (uv + cos_theta * n);
	// calc R`|| - fabs return absolut value of x, - to + 00
	Vec3f r_out_parallel = -sqrt(fabs(1.0f - r_out_perp.norm())) * n;

	return r_out_perp + r_out_parallel;

}

class material {
public:
	virtual bool scatter(const Ray& r_in, const hit_record& rec,
		Colour& attenuation, Ray& scattered)const = 0;
	virtual Colour emitted() const {
		return Colour(0, 0, 0);
	};
};

//class diffuse_light : public material {
//public:
//	diffuse_light() {}
//	diffuse_light(Colour c) : emit(make_shared<Colour>(c)) {}
//	virtual bool scatter(const Ray& r_in, const hit_record& rec,
//		Colour& attenuation, Ray& scattered)const override {
//		return false;
//	}
//
//	shared_ptr<Colour> emit;
//};

class matte : public material {
public:
	matte(const Colour& a) : albedo(a) {};

	virtual bool scatter(const Ray& r_in, const hit_record& rec,
		Colour& attenuation, Ray& scattered) const override {
		auto scatter_direction = -rec.normal ;

		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;
		// set incoming scattered and attenuation to the sphere albedo 
		// which is called in ray_color - another ray will be sent out till this will return false
		scattered = Ray(rec.p, scatter_direction);
		attenuation = albedo;
		return true;
	}

	Colour albedo;
};

class lambertian : public material {
public: 
	lambertian(const Colour& a) : albedo(a) {};

	virtual bool scatter(const Ray& r_in, const hit_record& rec, 
		Colour& attenuation, Ray& scattered ) const override{
		auto scatter_direction = rec.normal + Vec3f().random_in_unit_sphere();

		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;
		// set incoming scattered and attenuation to the sphere albedo 
		// which is called in ray_color - another ray will be sent out till this will return false
		scattered = Ray(rec.p, scatter_direction);
		attenuation = albedo;
		return true;
	}

	Colour albedo;
};

class lambertianT : public material {
public:

	lambertianT(const Colour& a) : albedo(make_shared<solid_color>(a)) {};
	lambertianT(shared_ptr<texture> a) : albedo(a) {};

	virtual bool scatter(const Ray& r_in, const hit_record& rec,
		Colour& attenuation, Ray& scattered) const override {
		auto scatter_direction = rec.normal + Vec3f().random_in_unit_sphere(); //random_unit_vector()

		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;
		// set incoming scattered and attenuation to the sphere albedo 
		// which is called in ray_color - another ray will be sent out till this will return false
		scattered = Ray(rec.p, scatter_direction);
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}

	shared_ptr<texture> albedo;
};


class metal : public material {
public: 
	metal(const Colour& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
	
	virtual bool scatter(const Ray& r_in, const hit_record& rec,
		Colour& attenuation, Ray& scattered) const override {

		Vec3f reflected = reflect(r_in.direction().normalize(), rec.normal);
		// shoot another from point of intersaction in the reflection vector direction + some fuzzines if there is any
		scattered = Ray(rec.p, reflected+fuzz*Vec3f().random_in_unit_sphere());
		attenuation = albedo;

		// calculate it`s direction 
		return (scattered.direction().dotProduct(rec.normal) > 0);
	}

	Colour albedo;
	double fuzz;
};

class diaelectric : public material {
public:
	diaelectric(double index_of_refraction) : ir(index_of_refraction){}

	virtual bool scatter(const Ray& r_in, const hit_record& rec,
		Colour& attenuation, Ray& scattered) const override {

		attenuation = Colour(1.0, 1.0, 1.0f);
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

		// vector needs to be normalized for dot prod to work and get cos0 in refract 
		Vec3f unit_direction = r_in.direction().normalize();

		double cos_theta = fmin(-unit_direction.dotProduct(rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

		// if the ray is inside the glass q = 1.5 and (air)q` = 1; sin0 1.5/1*sin0 >1.0 
		// case where glass cannot refract and will be reflected
		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		Vec3f direction;

		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract(unit_direction, rec.normal, refraction_ratio);

		scattered = Ray(rec.p, direction);
		
		return true;
	}

	double ir;

private:
	static double reflectance(double cosine, double ref_idx) {
		// Shlick`s approximation for reflectance
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}
};
//
//class metal : public material {
//public:
//	metal(const Colour& a) : albedo(a) {}
//
//	virtual bool scatter(const Ray& r_in, const hit_record& rec,
//		Colour& attenuation, Ray& scattered) const override {
//
//		Vec3f reflected = Vec3f::reflect(r_in.direction().normalize(), rec.normal);
//
//		scattered = Ray(rec.p, reflected );
//		attenuation = albedo;
//		return (scattered.direction().dotProduct(rec.normal) > 0);
//	}
//
//	Colour albedo;
//	double fuzz;
//};
//

class diffuse_light : public material {
public:
	diffuse_light() {}
	diffuse_light(Colour c) : emit(make_shared<Colour>(c)) {};

	virtual bool scatter(const Ray& r_in, const hit_record& rec,
		Colour& attenuation, Ray& scattered) const override {
		return false;
	}

	virtual Colour emitted() const override {
		return *emit;
	}

	shared_ptr<Colour> emit;
};

