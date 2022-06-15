#pragma once

#include "geometry.h"
#include "tgaimage.h"
#include <vector>
#include <SDL.h>

class Model {
private:
	std::vector<Vec3f> verts_;              // Stores Vec3f for every model vertex world position
	std::vector<std::vector<int> > faces_;
	std::vector<std::vector<int>> vnorm_; // Stores a vector of vector<int> that represent indices in verts_ for vertices comprising a face
	std::vector<Vec2f> vts_;				// Stores Vec3f for every model vertex texture coordinate
	std::vector<Vec3f> vns_;

public:
	Model(const char *filename);
	~Model();
	int nverts();
	int nfaces();
	Vec3f vert(int i);
	Vec3f vn(int i);
	std::vector<int> vnorm(int i);
	Vec2f vt(int i);
	std::vector<int> face(int idx);
	std::string name;
	TGAImage texture;

};

