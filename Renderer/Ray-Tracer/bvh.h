#pragma once
#include "common.h"
#include "Hittable.h"
#include "Hittable_list.h"

class bvh_node : public Hittable
{
public:
	bvh_node() {};

	bvh_node(const Hittable_list& list): bvh_node(list.objects, 0 , list.objects.size()){}
	
	// container of boxes
	bvh_node(const std::vector < shared_ptr < Hittable>>& src_objects, size_t start,
		size_t end);
	
	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

	virtual bool bounding_box(AABB& output_box) const override;

public:
	// left and right pointers to generic hittables (any primitive) allow us to split hierarcy
	shared_ptr<Hittable> left;
	shared_ptr<Hittable> right;
	AABB box;
};

bool bvh_node::bounding_box(AABB& output_box) const {
	output_box = box;
	return true;
}

// does not check for null pointers (!container planes)
bool bvh_node::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
	// if container not hit
	if (!box.hit(r, t_min, t_max)) return false;

	bool hit_left = left->hit(r, t_min, t_max, rec);
	bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

	// if left is hit then checks the right else false 
	return hit_left || hit_right; 
}




inline int random_int(int min, int max) {
	return static_cast<int>(random_double(min, max + 1));
}

// return true if min < max bb point 
inline bool box_compare(const shared_ptr < Hittable> a, const shared_ptr<Hittable> b, int axis) {
	AABB box_a;
	AABB box_b;

	if (!a->bounding_box(box_a) || !b->bounding_box(box_b))
		std::cerr << "No bounding box in bvh_node constructor.\n";

	return box_a.min()[axis] < box_b.min()[axis];
}

bool box_x_compare(const shared_ptr<Hittable> a, const shared_ptr<Hittable> b) {
	return box_compare(a, b, 0);
}

bool box_y_compare(const shared_ptr<Hittable> a, const shared_ptr<Hittable> b) {
	return box_compare(a, b, 1);
}

bool box_z_compare(const shared_ptr<Hittable> a, const shared_ptr<Hittable> b) {
	return box_compare(a, b, 2);
}

bvh_node::bvh_node(const std::vector < shared_ptr < Hittable>>& src_objects, size_t start,
	size_t end) {
	// all list of hittable
	auto objects = src_objects;

	int axis = random_int(0, 2);

	auto comparator = (axis == 0) ? box_x_compare :
		(axis == 1) ? box_y_compare :
		box_z_compare;

	// constructor: end - size of list - default start = 0
	size_t object_span = end - start;

	if (object_span == 1) {
		left = right = objects[start];
	}
	// small bb point = left 
	else if (object_span == 2) {
		if (comparator(objects[start], objects[start + 1])) {
			left = objects[start];
			right = objects[start + 1];
		}
		else
		{
			left = objects[start + 1];
			right = objects[start];
		}
	}
	else
	{
		std::sort(objects.begin() + start, objects.begin() + end, comparator);
		auto mid = start + object_span / 2;
		left = make_shared<bvh_node>(objects, start, mid);
		right = make_shared<bvh_node>(objects, mid, end);
	}

	AABB box_left, box_right;

	if (!left->bounding_box(box_left) || !right->bounding_box(box_right))
		std::cerr << "No bounding box in bvh_node constructor.\n";

	box = surrounding_box(box_left, box_right);
}
