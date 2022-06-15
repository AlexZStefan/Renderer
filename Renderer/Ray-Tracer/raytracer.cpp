#pragma once

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "SDL.h" 
#include "Hittable.h"
#include "Sphere.h"
#include "Hittable_list.h"
#include "common.h"
#include "Ray.h"
#include "Camera.h"
#include "ThreadPool.h"
#include "TriangleIntersection.h"
#include "bvh.h"

#include <fstream>
#include <chrono>

#define M_PI 3.14159265359

SDL_Window* window;
SDL_Renderer* renderer;
// 480 by 640
SDL_Surface* screen;
TGAImage scene_image; 
TGAColor scene_color;
int const imageWidth = 1920;
int const imageHeight = 1080;

void init() {
	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window* window = SDL_CreateWindow(
		"Software Ray Tracer",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		1920,
		1080,
		0
	);

	screen = SDL_GetWindowSurface(window);

	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
}

void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
	int bpp = surface->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to set */
	Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

	switch (bpp) {
	case 1:
		*p = pixel;
		break;

	case 2:
		*(Uint16*)p = pixel;
		break;

	case 3:
		if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
			p[0] = (pixel >> 16) & 0xff;
			p[1] = (pixel >> 8) & 0xff;
			p[2] = pixel & 0xff;
		}
		else {
			p[0] = pixel & 0xff;
			p[1] = (pixel >> 8) & 0xff;
			p[2] = (pixel >> 16) & 0xff;
		}
		break;

	case 4:
		*(Uint32*)p = pixel;
		break;
	}
}

// method to ensure colours don’t go out of 8 bit range in RGB​
void clamp255(Vec3f& col) {
	col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
	col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
	col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}


// mages color, reflectance, refractance of the ray - based on material
Colour ray_colour(const Ray& r, const Colour & background, const Hittable& world, int depth) {
	hit_record rec;

	// if exeded ray bounce limit, no more light is gathered
	if (depth <= 0) return Colour(0, 0, 0);

	// check all objects in the world against the ray and if intersected then returns the color gatherd from the ray 
	// else returns background 
	if (!world.hit(r, 0.001, infinity, rec))return background;
	
		Ray scattered;
		Colour attenuation;
		Colour emitted = rec.mat_ptr->emitted();

		//  called till return false - attenuation = material albedo or 1,1,1 if diaelectric 
		// scatter checks material and reflects/refracts - passed by refference  
		// checked till depth is 0 - if so return Color: 0 0 0 
		// else shoot another ray in the world and acumulate color
		if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return emitted;

		return emitted + attenuation * ray_colour(scattered,background, world, depth - 1);

}


auto mat_glass = make_shared<diaelectric>(1.4);
auto mat_table = make_shared<lambertian>(Colour(1, 1, 1));
auto mat_wall = make_shared<matte>(Colour(.90, .90, .90));
auto mat_tray = make_shared<lambertian>(Colour(0.4, .2, .2));
auto mat_plate = make_shared<lambertian>(Colour(0.9, .9, .9));
auto mat_coffaccess = make_shared<lambertian>(Colour(0, 0, 0));
auto mat_coffbody = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0);
auto mat_metal = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0);
auto mat_mirror = make_shared<metal>(Colour(0.9, 0.9, 0.9), 0);
auto mat_mug = make_shared<metal>(Colour(0.4, 0.7, 0.5), 1);
auto mat_mirror_frame = make_shared<metal>(Colour(.7, 0.3, 0.2), 1);
auto mat_spoon = make_shared<metal>(Colour(0.5, 0.5, 0.5), 0);
auto mat_trayr = make_shared<metal>(Colour(1.0, 0.6, 0.6), 1);

auto orangeText = make_shared<image_texture>("orange_text.tga");
auto mat_diffuse_light = make_shared<diffuse_light>(Colour(255,255,255));

Hittable_list random_scene() {
	Hittable_list world;

	auto ground_material = make_shared<lambertian>(Colour(0.5, 0.5, 0.5));
	auto checker = make_shared<checker_texture>(Colour(0.2, 0.3, 0.1), Colour(0.9, 0.9, 0.9));

	auto earth_texture = make_shared<image_texture>("earthmap.tga");
	auto earth_surface = make_shared<lambertianT>(earth_texture);
	world.add(make_shared<Sphere>(Point3f(0, -100, 0), 100, earth_surface));

	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			auto chose_mat = random_double();
			Point3f centre(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

			if ((centre - Point3f(4, 0.2, 0)).length() > 0.9) {
				shared_ptr<material> sphere_material;

				if (chose_mat < 0.8) {
					auto albedo = Colour::random() * Colour::random();
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
				}
				else if(chose_mat< 0.95)
				{
					auto albedo = Colour::random(0.5, 1);
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
				}
				else {
					sphere_material = make_shared<diaelectric>(1.5);
					world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
				}
			}
		}
	}

	auto material = make_shared<diaelectric>(1.5);
	auto material1 = make_shared<lambertian>(Colour(0.4, 0.2, 0.1));
	auto material2 = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.0);

	world.add(make_shared<Sphere>(Point3f(0, 1, 0), 1.0, material));
	world.add(make_shared<Sphere>(Point3f(-4, 1, 0), 1.0, material1));
	world.add(make_shared<Sphere>(Point3f(4, 1, 0), 1.0, material2));

	//return world;
	return Hittable_list(make_shared<bvh_node>(world));
	return world;

}

Hittable_list test_scene(std::vector<Model*>models) {
	Hittable_list world;

	Vec3f transform = (0, 0, 0);

	for (auto& m : models)
	{
		for (uint32_t i = 0; i < m->nfaces(); ++i) {
			const Vec3f& v0 = m->vert(m->face(i)[0]);
			const Vec3f& v1 = m->vert(m->face(i)[1]);
			const Vec3f& v2 = m->vert(m->face(i)[2]);

			const Vec3f& vn0 =m->vn(m->vnorm(i)[0]);
			const Vec3f& vn1 =m->vn(m->vnorm(i)[1]);
			const Vec3f& vn2 =m->vn(m->vnorm(i)[2]);
			if (m->name == "mplate")
			{
				world.add(make_shared<TriangleIntersection>(v0 , v1 , v2 , vn0, vn1, vn2, mat_plate));
			}
			else if (m->name == "mtable")
			{
				world.add(make_shared<TriangleIntersection>(v0 , v1 , v2, vn0, vn1, vn2, mat_table));
			}
			else if (m->name == "mmug")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_mug));
			}
			else if (m->name == "mmirror")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_mirror));
			}
			
			else if (m->name == "mframemirror")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_mirror_frame));
			}
			else if (m->name == "mcoffbody")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_coffbody));
			}
			else if (m->name == "mspoon")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_spoon));
			}
			else if (m->name == "mtray")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_tray));
			}
			else if (m->name == "mglass")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_glass));
			}
			else if (m->name == "mwall")
			{
				world.add(make_shared<TriangleIntersection>(v0, v1, v2, vn0, vn1, vn2, mat_wall));
			}
			
			else
			{
				world.add(make_shared<TriangleIntersection>(v0, v1 , v2 , vn0, vn1, vn2, mat_coffaccess));

			}
		}
	}

	return Hittable_list(make_shared<bvh_node>(world));
}


// handles shooting the ray from camera - through lens to all objects in hittable list 
// set pixels on the screen after max_depth reaches  0
// y - hittable_list, pix-colour, camera, nrOfSamples, image_width, image_height, depth, scale
void LineRenderer(int &y, Hittable_list &world, Camera& cam, int& spp,
	const int& image_width, const int& image_height, const int &max_depth, float &scale) {
	Colour pix_col;
	Colour background;

	for (int x = 0; x < screen->w; ++x) {
		pix_col = Colour(0, 0, 0); // reset col per pix
		// samples each pixel 
		for (int s = 0; s < spp; s++) {
			auto u = double(x + random_double()) / (image_width - 1);
			auto v = double(y + random_double()) / (image_height - 1);
			// new ray with cam.origin u and v direction 
			Ray ray = cam.get_ray(u, v);
			Vec3f unit_direction = ray.direction().normalize();
			auto t = 0.5 * (unit_direction.y + 1.0);
			background = (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(1,1, 1.0) * 255;//0.5,0.7,1
			// compute the pixel color from the ray shot in the world 
			pix_col = pix_col + ray_colour(ray, background, world, max_depth);

		}
		pix_col = (pix_col / 255.f) * spp;
		pix_col.x = sqrt(pix_col.x);
		pix_col.y = sqrt(pix_col.y);
		pix_col.z = sqrt(pix_col.z);
		pix_col *= 255;

		Uint64 colour = SDL_MapRGB(screen->format, pix_col.x * scale, pix_col.y * scale, pix_col.z * scale);
		putpixel(screen, x, y, colour);
	}
}

// SDL get pixel on surface 
Uint32 getpixel(SDL_Surface* surface, int x, int y)
{
	//source from Danny (https://stackoverflow.com/questions/53033971/how-to-get-the-color-of-a-specific-pixel-from-sdl-surface) starts here
	int bpp = surface->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to retrieve */
	Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

	switch (bpp)
	{
	case 1:
		return *p;
		break;

	case 2:
		return *(Uint16*)p;
		break;

	case 3:
		if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
			return p[0] << 16 | p[1] << 8 | p[2];
		else
			return p[0] | p[1] << 8 | p[2] << 16;
		break;

	case 4:
		return *(Uint32*)p;
		break;

	default:
		return 0;       /* shouldn't happen, but avoids warnings */
	}


	/* source from Danny ends here */

}

int main(int argc, char** argv)
{
	// initialise SDL2
	init();
	// create pool from number of existing threads. 

	std::vector<Model*> models; 

	Model* mtable = new Model("resources/mtable.obj");
	Model* mplate = new Model("resources/mplate.obj");
	Model* mcoffhand = new Model("resources/mcoffhand.obj");
	Model* mcoffbody = new Model("resources/mcoffbody.obj");
	Model* mcoffetop = new Model("resources/mcoffetop.obj");
	Model* mtray = new Model("resources/mtray.obj");
	Model* mwall = new Model("resources/mwall.obj");
	Model* mspoon = new Model("resources/mspoon.obj");
	Model* mmirror = new Model("resources/mmirror.obj");
	Model* mframemirror = new Model("resources/mframemirror.obj");
	Model* mmug = new Model("resources/mmug.obj");
	Model* mglass = new Model("resources/mglass.obj");

	mmug->name = "mmug";
	mtable->name = "mtable";
	mplate->name = "mplate";
	mmirror->name = "mmirror";
	mframemirror->name = "mframemirror";
	mcoffbody->name = "mcoffbody";
	mspoon->name = "mspoon";
	mtray->name = "mtray";
	mglass->name = "mglass";
	mwall->name = "mwall";


	models.push_back(mspoon);
	models.push_back(mtable);
	models.push_back(mglass);
	models.push_back(mplate);
	models.push_back(mtray);
	models.push_back(mcoffetop);
	models.push_back(mcoffbody);
	models.push_back(mcoffhand);
	models.push_back(mmug);
	models.push_back(mmirror);
	models.push_back(mframemirror);
	models.push_back(mwall);
	

	int spp = 100;
	float scale = 1.0f / spp;
	const int max_depth = 5;
	//auto R = cos(pi / 4);

	Point3f lookFrom(0, 4, -8.6);
	Point3f lookAt(0, 0, 15);
	Vec3f vup(0, 1, 0);
	auto dist_to_focus = 8;//(lookFrom - lookAt).length();
	auto aperture = 0.05f;
	//16/9 // 0.61
	const auto aspect_ratio = 16 / 9;
	const int& image_width = screen->w;
	const int image_height = static_cast<int> (image_width / aspect_ratio);

	Camera cam(lookFrom, lookAt, vup, 60, aspect_ratio, aperture, dist_to_focus);
	Colour pix_col;

	scene_image = TGAImage(imageWidth, imageHeight, 3);

	//Scene
	Hittable_list world = test_scene(models);
	//Hittable_list world = random_scene();
	world.add(make_shared<Sphere>(Point3f(0, 20, 0), 1.0, mat_diffuse_light));

	auto orangeT = make_shared<lambertianT>(orangeText);
	world.add(make_shared<Sphere>(Point3f(0.2, 0.4, 3), 0.2, orangeT));
	world.add(make_shared<Sphere>(Point3f(-0.25, 0.4, 3.1), 0.2, orangeT));
	world.add(make_shared<Sphere>(Point3f(0.4, 0.4, 2.8), 0.2, orangeT));

	double t;
	SDL_Event e;
	bool running = true;
	auto t_start = std::chrono::high_resolution_clock::now();


		t_start = std::chrono::high_resolution_clock::now();

		//while (running)
		//{


		// clear back buffer, pixel data on surface and depth buffer (as movement)
		SDL_FillRect(screen, nullptr, SDL_MapRGB(screen->format, 0, 0, 0));
		SDL_RenderClear(renderer);

		
		// ray loop 
		{
			ThreadPool pool(std::thread::hardware_concurrency());

			for (int y = 0; y < screen->h - 1; y++) {
				//std::cerr << "\rScanlines remaining: " << y << std::flush;

				pool.enqueue(std::bind(&LineRenderer, y, world, cam, spp,
					image_width, image_height, max_depth, scale));
		
			}
			//std::this_thread::sleep_for(std::chrono::milliseconds(2000));
		}
		
		auto t_end = std::chrono::high_resolution_clock::now();
		auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
		std::cerr << "Frame render time:  " << passedTime << " ms" << std::endl;
		
		SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, screen);
		if (texture == NULL) {
			fprintf(stderr, "CreateTextureFromSurface failed: %s\n", SDL_GetError());
			exit(1);
		}

		SDL_RenderCopyEx(renderer, texture, nullptr, nullptr, 0, 0, SDL_FLIP_VERTICAL);
		SDL_RenderPresent(renderer);

		/*SDL_FreeSurface(screen);
		SDL_DestroyTexture(texture);			
		}*/

		
		int bpp = screen->format->BytesPerPixel;

		for (int x = 0; x < screen->w; x++)
		{
			for (int y = 0; y < screen->h; y++) {

				SDL_Color rgb;
				Uint32 data = getpixel(screen, x, y);
				SDL_GetRGB(data, screen->format, &rgb.r, &rgb.g, &rgb.b);
				scene_color = TGAColor(rgb.r, rgb.g, rgb.b, bpp);
				scene_image.set(x, y, scene_color);
			}
		}

		scene_image.write_tga_file("output.tga");
		SDL_FreeSurface(screen);
		SDL_DestroyTexture(texture);

		if (SDL_PollEvent(&e))
		{
			switch (e.type) {
			case SDL_KEYDOWN:
				switch (e.key.keysym.sym) {
				case SDLK_ESCAPE:

					int bpp = screen->format->BytesPerPixel;

					for (int x = 0; x < imageWidth; x++)
					{
						for (int y = 0; y < imageWidth; y++) {

							SDL_Color rgb;
							Uint32 data = getpixel(screen, x, y);
							SDL_GetRGB(data, screen->format, &rgb.r, &rgb.g, &rgb.b);
							scene_color = TGAColor(rgb.r, rgb.g, rgb.b, bpp);
							scene_image.set(x, y, scene_color);
						}
					}

					scene_image.write_tga_file("output.tga");

					


					running = false;
					break;
				}
				break;
			}
		}
		for (int i = models.size() - 1; i > -1; i--) {
			delete models[i];
		}

	
	

	return 0;
}
