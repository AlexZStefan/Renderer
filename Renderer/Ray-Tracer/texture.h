#pragma once
#include "geometry.h"

class texture
{
public:
	virtual Colour value(double u, double v, const Point3f& p) const = 0;
};

class solid_color : public texture {

public:
	solid_color(){}
	solid_color(Colour c): color_value(c){}

	solid_color(double red, double green, double blue): solid_color(Colour(red,green,blue)){}

	virtual Colour value(double u, double v, const Point3f& p) const override {
		return color_value;
	}

private:
	Colour color_value;

};

class checker_texture : public texture {
public:
	checker_texture(){}
	checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd) : even(_even), odd(_odd){}

	checker_texture(Colour c1, Colour c2) : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {};

	virtual Colour value(double u, double v, const Point3f& p) const override {
		auto sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);

		if (sines < 0)
			return odd->value(u, v, p);
		else
		{
			return even->value(u, v, p);
		}
	}

	shared_ptr<texture> odd;
	shared_ptr<texture> even;
};

class image_texture : public texture {
public:

	image_texture() :  width(0), height(0), bytes_per_scanline(0){}

	image_texture(const char* filename) {
		img.read_tga_file(filename);
		int bytes_per_pixel = img.get_bytespp();
		auto components_per_pixel = bytes_per_pixel;
		width = img.get_width();
		height = img.get_height();
		components_per_pixel = img.get_bytespp();

		bytes_per_scanline = components_per_pixel * width;
	}

	~image_texture() {
	}

	// for pixel(i, j) in an Nx by Ny image, the image texture position is 
	// u = i/ Nx -1 ; v = j/Ny-1
	virtual Colour value(double u, double v, const Vec3f& p) const override {
		// return cyan for debugging purps
		if (img.data == nullptr) return Colour(0, 1, 1);

		u = clamp(u, 0, 1);
		v = 1.0 - clamp(v, 0, 1); // flip v to img coord 

		auto i = static_cast<int>(u * width);
		auto j = static_cast<int>(v * height);

		if (i >= width) i = width - 1;
		if (j >= height) i = height - 1;

		const auto color_scale = 1.0 / 255.0;
		auto pixel = img.data + j * bytes_per_scanline + i * img.bytespp;

		TGAColor imgC = TGAColor(img.data + j * bytes_per_scanline + i * img.bytespp, img.bytespp);
		//auto d = data[1][7];
		//auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

		return Colour(color_scale * imgC.r, color_scale * imgC.g, color_scale * imgC.b);
	}

private:
	int width, height; 
	int bytes_per_scanline; 
	TGAImage img;
};