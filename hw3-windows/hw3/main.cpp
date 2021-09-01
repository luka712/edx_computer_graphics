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
#include <sstream>
#include <stack>
#include "include/FreeImage.h"
#include "bones_lights.hpp"



// TODO: move to raycaster 
struct Scene
{
	std::vector<bns::Camera> Cameras;
	std::vector<bns::BaseShape*> Shapes;
	std::vector<bns::BaseLight*> Lights;
};

void CreateSceneFromFile(const char* filename, Scene& scene)
{
	std::vector<bns::Vec3F> vertices;

	I32 camera_width = 0;
	I32 camera_height = 0;

	bns::FileContents* file = bns::ReadAndCloseFile(filename);

	std::string line;
	std::stringstream file_stream(file->Contents);


	std::stack<bns::Mat4x4F> matrix_stack;
	matrix_stack.push(bns::Mat4x4F::Identity());

	bns::ColorF ambient_color = bns::ColorF(.2f, .2f, .2f);
	bns::ColorF emission_color;
	bns::ColorF diffuse_color;
	bns::ColorF specular_color;
	F32 shininess = 0.0f;
	F32 attenuation_constant = 1.0f;
	F32 attenuation_linear = 0.0f;
	F32 attenuation_quadratic = 0.0f;

	while (std::getline(file_stream, line))
	{
		if (line.rfind("end", 0) == 0)
		{
			break;
		}
		if (line.rfind("size", 0) == 0)
		{
			I32 end_index = 0;

			ReadI32FromString(line.c_str(), end_index, &end_index, &camera_width);
			ReadI32FromString(line.c_str(), end_index, &end_index, &camera_height);
		}
		else if (line.rfind("camera", 0) == 0)
		{
			I32 end_index = 0;
			const char* ptr_to_camera = line.c_str();

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

			scene.Cameras.push_back(bns::Camera(eyeinit, center, upinit, fov, camera_width, camera_height));
		}
		else if (line.rfind("vertex", 0) == 0)
		{
			const char* ptr_to_vertex = line.c_str();

			I32 end_index = 0;
			bns::Vec3F vert;
			ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.X);
			ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.Y);
			ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.Z);
			vertices.push_back(vert);

		}
		else if (line.rfind("pushTransform", 0) == 0)
		{
			if (matrix_stack.empty()) 
			{
				matrix_stack.push(bns::Mat4x4F::Identity());
			}
			else
			{
				bns::Mat4x4F m = matrix_stack.top();
				m =  bns::Mat4x4F::Identity() * m;
				matrix_stack.push(m);
			}
		}
		else if (line.rfind("popTransform", 0) == 0)
		{
			matrix_stack.pop();
		}
		else if (line.rfind("translate", 0) == 0)
		{
			const char* ptr = line.c_str();

			I32 end_index = 0;
			bns::Vec3F v;
			ReadF32FromString(ptr, end_index, &end_index, &v.X);
			ReadF32FromString(ptr, end_index, &end_index, &v.Y);
			ReadF32FromString(ptr, end_index, &end_index, &v.Z);

			bns::Mat4x4F m = matrix_stack.top();
			matrix_stack.pop();
			m = m * bns::Mat4x4F::Translate(v);
			matrix_stack.push(m);
		}
		else if (line.rfind("scale", 0) == 0)
		{
			const char* ptr = line.c_str();

			I32 end_index = 0;
			bns::Vec3F v;
			ReadF32FromString(ptr, end_index, &end_index, &v.X);
			ReadF32FromString(ptr, end_index, &end_index, &v.Y);
			ReadF32FromString(ptr, end_index, &end_index, &v.Z);

			bns::Mat4x4F m = matrix_stack.top();
			matrix_stack.pop();
			m = m * bns::Mat4x4F::Scale(v);
			matrix_stack.push(m);
		}
		else if (line.rfind("rotate", 0) == 0)
		{
			const char* ptr = line.c_str();

			I32 end_index = 0;
			bns::Vec3F v;
			F32 theta;
			ReadF32FromString(ptr, end_index, &end_index, &v.X);
			ReadF32FromString(ptr, end_index, &end_index, &v.Y);
			ReadF32FromString(ptr, end_index, &end_index, &v.Z);
			ReadF32FromString(ptr, end_index, &end_index, &theta);

			bns::Mat4x4F m = matrix_stack.top();
			matrix_stack.pop();
			m = m * bns::Mat4x4F::RotationMatrix(bns::Radians(theta), v);
			matrix_stack.push(m);
		}
		else if (line.rfind("ambient", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			ReadF32FromString(ptr, end_index, &end_index, &ambient_color.R);
			ReadF32FromString(ptr, end_index, &end_index, &ambient_color.G);
			ReadF32FromString(ptr, end_index, &end_index, &ambient_color.B);

		}
		else if (line.rfind("emission", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			ReadF32FromString(ptr, end_index, &end_index, &emission_color.R);
			ReadF32FromString(ptr, end_index, &end_index, &emission_color.G);
			ReadF32FromString(ptr, end_index, &end_index, &emission_color.B);
		}
		else if (line.rfind("diffuse", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			ReadF32FromString(ptr, end_index, &end_index, &diffuse_color.R);
			ReadF32FromString(ptr, end_index, &end_index, &diffuse_color.G);
			ReadF32FromString(ptr, end_index, &end_index, &diffuse_color.B);
		}
		else if (line.rfind("specular", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			ReadF32FromString(ptr, end_index, &end_index, &specular_color.R);
			ReadF32FromString(ptr, end_index, &end_index, &specular_color.G);
			ReadF32FromString(ptr, end_index, &end_index, &specular_color.B);
		}
		else if (line.rfind("shininess", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			ReadF32FromString(ptr, end_index, &end_index, &shininess);
		}
		else if (line.rfind("attenuation", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			ReadF32FromString(ptr, end_index, &end_index, &attenuation_constant);
			ReadF32FromString(ptr, end_index, &end_index, &attenuation_linear);
			ReadF32FromString(ptr, end_index, &end_index, &attenuation_quadratic);


		}
		else if (line.rfind("tri", 0) == 0)
		{
			const char* ptr = line.c_str();

			I32 end_index = 0;
			I32 i_1 = 0;
			I32 i_2 = 0;
			I32 i_3 = 0;
			bns::Vec3F vert;
			ReadI32FromString(ptr, end_index, &end_index, &i_1);
			ReadI32FromString(ptr, end_index, &end_index, &i_2);
			ReadI32FromString(ptr, end_index, &end_index, &i_3);

			bns::Vec3F A = vertices[i_1];
			bns::Vec3F B = vertices[i_2];
			bns::Vec3F C = vertices[i_3];

			scene.Shapes.emplace_back(new bns::TriangleShape());
			bns::TriangleShape* ptr_triangle = static_cast<bns::TriangleShape*>(scene.Shapes.back());
			ptr_triangle->Triangle = { A, B, C };
			ptr_triangle->SetTransform(matrix_stack.top());
			ptr_triangle->Material.Ambient = ambient_color;
			ptr_triangle->Material.Emission = emission_color;
			ptr_triangle->Material.Diffuse = diffuse_color;
			ptr_triangle->Material.Specular = specular_color;
			ptr_triangle->Material.Shininess = shininess;
		}
		else if (line.rfind("sphere", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;
			F32 x = 0;
			F32 y = 0;
			F32 z = 0;
			F32 r = 0; // radius
			bns::Vec3F vert;
			ReadF32FromString(ptr, end_index, &end_index, &x);
			ReadF32FromString(ptr, end_index, &end_index, &y);
			ReadF32FromString(ptr, end_index, &end_index, &z);
			ReadF32FromString(ptr, end_index, &end_index, &r);

			scene.Shapes.emplace_back(new bns::SphereShape());
			bns::SphereShape* ptr_sphere = static_cast<bns::SphereShape*>(scene.Shapes.back());
			ptr_sphere->Sphere.Position = bns::Vec3F(x, y, z);
			ptr_sphere->Sphere.Radius = r;
			ptr_sphere->SetTransform(matrix_stack.top());
			ptr_sphere->Material.Ambient = ambient_color;
			ptr_sphere->Material.Diffuse = diffuse_color;
			ptr_sphere->Material.Emission = emission_color;
			ptr_sphere->Material.Specular = specular_color;
			ptr_sphere->Material.Shininess = shininess;
			ptr_sphere->Material.RefractiveIndex = 1.5f;
		}
		else if (line.rfind("point", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			scene.Lights.emplace_back(new bns::PointLight());
			bns::PointLight& light = *(static_cast<bns::PointLight*>(scene.Lights.back()));
			ReadF32FromString(ptr, end_index, &end_index, &light.Position.X);
			ReadF32FromString(ptr, end_index, &end_index, &light.Position.Y);
			ReadF32FromString(ptr, end_index, &end_index, &light.Position.Z);
			ReadF32FromString(ptr, end_index, &end_index, &light.Color.R);
			ReadF32FromString(ptr, end_index, &end_index, &light.Color.G);
			ReadF32FromString(ptr, end_index, &end_index, &light.Color.B);

			light.Attenuation.Constant = attenuation_constant;
			light.Attenuation.Linear = attenuation_linear;
			light.Attenuation.Quadratic = attenuation_quadratic;
		}
		else if (line.rfind("directional", 0) == 0)
		{
			const char* ptr = line.c_str();
			I32 end_index = 0;

			scene.Lights.emplace_back(new bns::DirectionalLight());
			bns::DirectionalLight& light = *(static_cast<bns::DirectionalLight*>(scene.Lights.back()));
			ReadF32FromString(ptr, end_index, &end_index, &light.Direction.X);
			ReadF32FromString(ptr, end_index, &end_index, &light.Direction.Y);
			ReadF32FromString(ptr, end_index, &end_index, &light.Direction.Z);
			ReadF32FromString(ptr, end_index, &end_index, &light.Color.R);
			ReadF32FromString(ptr, end_index, &end_index, &light.Color.G);
			ReadF32FromString(ptr, end_index, &end_index, &light.Color.B);
		}
	}


	bns::FreeFileContents(file);
}

void saveScreenshot(std::string fname, void* pixels, I32 w, I32 h)
{
	FIBITMAP* img = FreeImage_ConvertFromRawBits((BYTE*)pixels, w, h, w * 4, 32, 0xFF000000, 0x00FF0000, 0x0000FF00, 0x000000FF);

	std::cout << "Saving screenshot: " << fname << "\n";

	FreeImage_Save(FIF_PNG, img, fname.c_str(), 0);
}



#define IMAGE 0

int main()
{
#if IMAGE

	FreeImage_Initialise(0);

	const char* source[2] = {
		//"scene4-ambient.test",
		//"scene4-diffuse.test",
		//"scene4-emission.test",
		//"scene4-specular.test",
		"scene5.test",
		"scene6.test",
		//"scene7.test"
	};
	const char* target[2] = {
		//"scene4-ambient.png",
		//"scene4-diffuse.png",
		//"scene4-emission.png",
		//"scene4-specular.png",
		"scene5.png",
		"scene6.png",
		//"scene7.png"
	};

	for (size_t i = 0; i < 1; i++)
	{
		Scene renderer_scene;
		CreateSceneFromFile(source[i], renderer_scene);

		I32 _index = 1;
		bns::Camera cam = renderer_scene.Cameras[0];


		bns::ColorF** colors = bns::AllocateColors(cam);
		bns::ThreadedRayTrace(cam,
			renderer_scene.Shapes.data(), renderer_scene.Shapes.size(),
			renderer_scene.Lights.data(), renderer_scene.Lights.size(),
			colors);
		void* pixels = bns::ColorsToARGB8888Pixels(cam, colors);


		saveScreenshot(std::string(target[i]), pixels, cam.ScreenWidth, cam.ScreenHeight);
		_index++;


		bns::FreeColors(cam, colors);
		bns::FreePixels(pixels);
	}

	FreeImage_DeInitialise();

#else 

	I32 _index = 1;

	Scene renderer_scene;
	CreateSceneFromFile("scene6.test", renderer_scene);

	bns::Camera cam = renderer_scene.Cameras[0];

	bns::ColorF** colors = bns::AllocateColors(cam);
	bns::ThreadedRayTrace(cam,
		renderer_scene.Shapes.data(), renderer_scene.Shapes.size(),
		renderer_scene.Lights.data(), renderer_scene.Lights.size(), colors);
	void* pixels = bns::ColorsToABGR8888Pixels(cam, colors);
	bns::FreeColors(cam, colors);

	for (auto it = renderer_scene.Shapes.begin(); it != renderer_scene.Shapes.end(); ++it)
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

	window = SDL_CreateWindow("Shooter 01", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, cam.ScreenWidth, cam.ScreenHeight, 0);

	if (!window)
	{
		printf("Failed to open %d x %d window: %s\n", cam.ScreenWidth, cam.ScreenHeight, SDL_GetError());
		exit(1);
	}

	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");

	renderer = SDL_CreateRenderer(window, -1, 0);

	if (!renderer)
	{
		printf("Failed to create renderer: %s\n", SDL_GetError());
		exit(1);
	}

	SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, cam.ScreenWidth, cam.ScreenHeight);
	SDL_Rect rect;
	rect.x = 0;
	rect.y = 0;
	rect.w = cam.ScreenWidth;
	rect.h = cam.ScreenHeight;
	SDL_UpdateTexture(texture, &rect, pixels, cam.ScreenWidth * 4);
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