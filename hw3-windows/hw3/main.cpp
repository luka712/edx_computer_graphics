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

struct Camera
{
	bns::Vec3F LookFrom;
	bns::Vec3F LookAt;
	bns::Vec3F up;
	F32 FOV;
};

template<typename T>
struct Intersection
{
	F32 MinDist;
	T* HitObject;
};

std::vector<bns::Vec3F> vertices;
std::vector<bns::TriangleF> triangles;

//bns::Vec3F upvector(const bns::Vec3F& up, const bns::Vec3F& zvec)
//{
//	bns::Vec3F x = bns::Vec3F::Cross(up, zvec);
//	bns::Vec3F y = bns::Vec3F::Cross(zvec, x);
//	y.Normalize();
//	return y;
//}

bns::RayF RayThroughPixel(Camera cam, I32 i, I32 j, I32 width, I32 height)
{
	// Steps
	// 1. Create a coordinate frame for the camera
	// 2. Define a rotation matrix
	// 3. Apply appropriate translation for camera ( eye ) location

	//      a          b x w
	// w = ---    u = -------       v = w x u
	//    ||a||     || b x w ||

	// a = eye - center
	bns::Vec3F eye = cam.LookFrom;
	bns::Vec3F center = cam.LookAt;
	bns::Vec3F up = cam.up;

	bns::Vec3F a = eye - center;
	bns::Vec3F w = bns::Vec3F::Normalize(a);

	bns::Vec3F b = bns::Vec3F::Normalize(up);

	bns::Vec3F b_cross_w = bns::Vec3F::Cross(b, w);
	bns::Vec3F b_cross_w_unit = bns::Vec3F::Normalize(b_cross_w);

	bns::Vec3F u = b_cross_w_unit;

	bns::Vec3F v = bns::Vec3F::Cross(w, u);

	F32 half_w = width / 2.0f;
	F32 half_h = height / 2.0f;

	F32 _i = i + 0.5f;
	F32 _j = j + 0.5f;

	F32 fov = bns::Radians(cam.FOV);
    
	// https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	F32 half_view = bns::Tan(fov / 2.0f);

	F32 aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
	
	// TODO: fix for aspect ration
	// in this assignemnt x is wider, therefore fix alpha, otherwise beta is to be fixed.

	F32 alpha = half_view * aspect_ratio * ((_i - half_w) / half_w);
	F32 beta = half_view *  ((half_h - _j) / half_h);

	// ray = eye + (alpha_u + beta_v - w)/||(alpha_u + beta_v - w)||

	bns::Vec3F dir = u * alpha + v * beta - w;

	dir.Normalize();

	bns::RayF ray(eye, dir);
	return ray;
}

Intersection<bns::TriangleF> GetIntersection(bns::RayF ray, std::vector<bns::TriangleF> triangles)
{
	F32 min_dist = MAX_F32;
	bns::TriangleF* ptr_hit_object = nullptr;
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		F32 t = ray.Intersects(*it);
		if (t > 0 && t < min_dist)
		{
			min_dist = t;
			ptr_hit_object = &*it;
		}
	}

	return {
		min_dist,
		ptr_hit_object
	};
}

bns::Vec3F FindColor(Intersection<bns::TriangleF> hit)
{
	// TODO: implement
	if (hit.MinDist != MAX_F32) 
	{
		return bns::Vec3F(1.0f, 0.0f, 0.0f);
	}
	return bns::Vec3F(0, 0, 0);
}

std::vector<std::vector<bns::Vec3F>> RayTrace(Camera cam, std::vector<bns::TriangleF> triangles, int width, int height)
{
	std::vector<std::vector<bns::Vec3F>> result(height);

	for (size_t j = 0; j < height; j++)
	{
		result[j] = std::vector<bns::Vec3F>(width);
		for (size_t i = 0; i < width; i++)
		{
			bns::RayF ray = RayThroughPixel(cam, i,j, width, height);
			Intersection<bns::TriangleF> hit = GetIntersection(ray, triangles);
			result[j][i] = FindColor(hit);
		}
	}

	return result;
}

int main()
{
	std::vector<Camera> cameras;

	I32 camera_width = 0;
	I32 camera_height = 0;

	bns::FileContents* file = bns::ReadAndCloseFile("scene1.test");

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

		cameras.push_back({ eyeinit, center, upinit,fov });

		ptr_to_camera = strstr(&ptr_to_camera[1], "camera");
	}

	char* ptr_to_maxverts = strstr(file->Contents, "maxverts");
	I32 maxverts = 0;
	ReadI32FromString(ptr_to_maxverts, 0, &maxverts);

	vertices = std::vector<bns::Vec3F>();


	char* ptr_to_vertex = strstr(file->Contents, "vertex");
	for (size_t i = 0; i < maxverts; i++)
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

		triangles.emplace_back(A, B, C);

		ptr_to_tri = strstr(&ptr_to_tri[1], "tri");
	}

	bns::FreeFileContents(file);

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
	

	/*
	SDL_Renderer* renderer;
	SDL_Window* window;

	int rendererFlags, windowFlags;

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

	SDL_Texture* texture =SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, camera_width, camera_height);
	SDL_Rect rect;
	rect.x = 0;
	rect.y = 0;
	rect.w = camera_width;
	rect.h = camera_height;
	SDL_UpdateTexture(texture,&rect,pixels,camera_width * 4);
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

	delete[] pixels;

	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_Quit();

	*/
	

	return 0;
}