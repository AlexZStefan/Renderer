#pragma once
#include "common.h"

class AABB
{
public:
	AABB(){}

	AABB(const Point3f& mini, const Point3f& maxi)
	{
		minimum = mini;
		maximum = maxi;
	}

	Point3f min()const { return minimum; }
	Point3f max()const { return maximum; }


	
	bool hit(const Ray& r, double t_min, double t_max) const {
		for (int i = 0; i < 3; i++) {
			auto t0 = fmin((minimum[i] - r.origin()[i]) / r.direction()[i], 
				(maximum[i] - r.origin()[i]) / r.direction()[i]);
			auto t1 = fmax((minimum[i] - r.origin()[i]) / r.direction()[i],
				(maximum[i] - r.origin()[i]) / r.direction()[i]);
			t_min = fmax(t0, t_min);
			t_max = fmin(t1, t_max);

			if (t_max <= t_min) return false;
		}
		return true;
	}

private:
	Point3f minimum;
	Point3f maximum;

	
};

// points that make up the surrounding box
AABB surrounding_box(AABB box0, AABB box1) {
	Point3f small(fmin(box0.min().x - 1e-3, box1.min().x - 1e-3),
		fmin(box0.min().y - 1e-3, box1.min().y - 1e-3),
		fmin(box0.min().z - 1e-3, box1.min().z - 1e-3));

	Point3f big(fmax(box0.max().x + 1e-3, box1.max().x + 1e-3),
		fmax(box0.max().y + 1e-3, box1.max().y + 1e-3),
		fmax(box0.max().z + 1e-3, box1.max().z + 1e-3));

	return AABB(small, big);
}