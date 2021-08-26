#include <cstdio>
#include <fstream>
#include <string>
#include "bones_types.hpp"
#include "bones_files.hpp"
#include <iostream>
#include "bones_char_string.hpp"
#include "bones_math.hpp"
#include <vector>
#include "loadbmp.h"
#define SDL_MAIN_HANDLED
#include <SDL.h>
#include "bones_raytracer.hpp"


std::vector<bns::Vec3F> vertices;
std::vector<bns::BaseShape*> shapes;
std::vector<bns::ColorF> ambient_colors;

U32 triangles_count = 0;
U32 spheres_count = 0;

int main()
{
	std::vector<bns::Camera> cameras;

	I32 camera_width = 0;
	I32 camera_height = 0;

	bns::FileContents* file = bns::ReadAndCloseFile("scene2.test");

	char* ptr_to_size = strstr(file->Contents, "size");
	if (ptr_to_size != 0)
	{
		I32 end_index = 0;
		// 5 to skip size + space
		ReadI32FromString(ptr_to_size, 4, &end_index, &camera_width);

		ReadI32FromString(ptr_to_size, end_index, &end_index, &camera_height);
	}

	char* ptr_to_camera = strstr(file->Contents, "camera");
	while (ptr_to_camera)
	{
		I32 end_index = 0;

		bns::Vec3F eyeinit;
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &eyeinit.X);
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &eyeinit.Y);
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &eyeinit.Z);

		bns::Vec3F center;
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &center.X);
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &center.Y);
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &center.Z);

		bns::Vec3F upinit;
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &upinit.X);
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &upinit.Y);
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &upinit.Z);

		F32 fov;
		ReadF32FromString(ptr_to_camera, end_index, &end_index, &fov);

		upinit = bns::Vec3F::Normalize(upinit);

		cameras.push_back(bns::Camera(eyeinit, center, upinit,fov,camera_width, camera_height)); 

		ptr_to_camera = strstr(&ptr_to_camera[1], "camera");
	}

	char* ptr_to_maxverts = strstr(file->Contents, "maxverts");
	I32 maxverts = 0;
	ReadI32FromString(ptr_to_maxverts, 0, &maxverts);

	vertices = std::vector<bns::Vec3F>();


	char* ptr_to_vertex = strstr(file->Contents, "vertex");
	for (I32 i = 0; i < maxverts; i++)
	{
		if (i > 0)
		{
			if (ptr_to_vertex != NULL)
			{
				ptr_to_vertex = strstr(&ptr_to_vertex[1], "vertex");
			}
		}
		I32 end_index = 0;
		bns::Vec3F vert;
		ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.X);
		ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.Y);
		ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.Z);
		vertices.push_back(vert);
	}


	char* ptr_to_tri = strstr(file->Contents, "tri");
	I32 index = 0;
	while (ptr_to_tri)
	{
		I32 end_index = 0;
		I32 i_1 = 0;
		I32 i_2 = 0;
		I32 i_3 = 0;
		bns::Vec3F vert;
		ReadI32FromString(ptr_to_tri, end_index, &end_index, &i_1);
		ReadI32FromString(ptr_to_tri, end_index, &end_index, &i_2);
		ReadI32FromString(ptr_to_tri, end_index, &end_index, &i_3);

		bns::Vec3F A = vertices[i_1];
		bns::Vec3F B = vertices[i_2];
		bns::Vec3F C = vertices[i_3];

		shapes.emplace_back(new bns::TriangleShape());
		bns::TriangleShape* ptr_triangle = static_cast<bns::TriangleShape*>(shapes.back());
		ptr_triangle->Triangle = { A, B, C };

		ptr_to_tri = strstr(&ptr_to_tri[1], "tri");
		triangles_count++;
	}

	char* ptr_to_ambient = strstr(file->Contents, "ambient");
	index = 0;
	while (ptr_to_ambient)
	{
		I32 end_index = 0;
		F32 r = 0;
		F32 g = 0;
		F32 b = 0;
		bns::Vec3F vert;
		ReadF32FromString(ptr_to_ambient, end_index, &end_index, &r);
		ReadF32FromString(ptr_to_ambient, end_index, &end_index, &g);
		ReadF32FromString(ptr_to_ambient, end_index, &end_index, &b);

		ambient_colors.emplace_back(r, g, b);

		ptr_to_ambient = strstr(&ptr_to_ambient[1], "ambient");
	}

	char* ptr_to_sphere = strstr(file->Contents, "sphere");
	index = 0;
	while (ptr_to_sphere)
	{
		I32 end_index = 0;
		F32 x = 0;
		F32 y = 0;
		F32 z = 0;
		F32 r = 0; // radius
		bns::Vec3F vert;
		ReadF32FromString(ptr_to_sphere, end_index, &end_index, &x);
		ReadF32FromString(ptr_to_sphere, end_index, &end_index, &y);
		ReadF32FromString(ptr_to_sphere, end_index, &end_index, &z);
		ReadF32FromString(ptr_to_sphere, end_index, &end_index, &r);

		shapes.emplace_back(new bns::SphereShape());
		bns::SphereShape* ptr_sphere = static_cast<bns::SphereShape*>(shapes.back());
		ptr_sphere->Sphere.Position = bns::Vec3F(x, y, z);
		ptr_sphere->Sphere.Radius = r;

		ptr_to_sphere = strstr(&ptr_to_sphere[1], "sphere");

		spheres_count++;
	}

	I32 ambient_index = 0;
	for (size_t i = 0; i < triangles_count; i+=2)
	{
		shapes[i]->Material.Color = ambient_colors[ambient_index];
		shapes[i+1]->Material.Color = ambient_colors[ambient_index++];
	}

	for (size_t i = triangles_count; i < triangles_count + spheres_count; i++)
	{
		shapes[i]->Material.Color = ambient_colors[ambient_colors.size() - 1];
	}

	bns::FreeFileContents(file);


#if IMAGE

	I32 _index = 1;
	for (auto it = cameras.begin(); it != cameras.end(); it++)
	{
		auto result = RayTrace(*it, triangles, camera_width, camera_height);
		unsigned char* pixels = (unsigned char*)malloc(camera_width * camera_height * 4);

		I32 pixel_index = 0;
		for (I32 y = 0; y < result.size(); y++)
		{
			for (I32 x = 0; x < result[y].size(); x++)
			{
				bns::Vec3F color = result[y][x];

				unsigned char r = static_cast<unsigned char>(color.X * 255);
				unsigned char g = static_cast<unsigned char>(color.Y * 255);
				unsigned char b = static_cast<unsigned char>(color.Z * 255);

				pixels[pixel_index] = r;
				pixels[pixel_index + 1] = g;
				pixels[pixel_index + 2] = b;
				pixels[pixel_index + 3] = 255;

				pixel_index += 4;
			}
		}

		unsigned int err = loadbmp_encode_file(_index + "_camera_test.bmp", pixels, camera_width, camera_height, LOADBMP_RGBA);
		_index++;

		free(pixels);
	}

#else 


	I32 _index = 1;

	bns::Camera cam = cameras[1];
	bns::ColorF** colors = bns::AllocateColors(cam);
	bns::RayTrace(cam, shapes.data(), shapes.size(), colors);
	void* pixels = bns::ColorsToABGR8888Pixels(cam, colors);
	bns::FreeColors(cam, colors);

	for (auto it = shapes.begin(); it != shapes.end(); ++it)
	{
		delete* it;
	}

	SDL_Renderer* renderer;
	SDL_Window* window;

	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) < 0)
	{
		printf("Couldn't initialize SDL: %s\n", SDL_GetError());
		exit(1);
	}

	window = SDL_CreateWindow("Shooter 01", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, camera_width, camera_height, 0);

	if (!window)
	{
		printf("Failed to open %d x %d window: %s\n", camera_width, camera_height, SDL_GetError());
		exit(1);
	}

	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");

	renderer = SDL_CreateRenderer(window, -1, 0);

	if (!renderer)
	{
		printf("Failed to create renderer: %s\n", SDL_GetError());
		exit(1);
	}

	SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, camera_width, camera_height);
	SDL_Rect rect;
	rect.x = 0;
	rect.y = 0;
	rect.w = camera_width;
	rect.h = camera_height;
	SDL_UpdateTexture(texture, &rect, pixels, camera_width * 4);
	SDL_RenderCopy(renderer, texture, &rect, &rect);
	SDL_RenderPresent(renderer);


	SDL_Event e;
	for (;;) {
		SDL_PollEvent(&e);
		if (e.type == SDL_QUIT) {
			SDL_Log("Program quit after %i ticks", e.quit.timestamp);
			break;
		}
	}

	bns::FreePixels(pixels);

	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_Quit();

#endif 

	return 0;
}