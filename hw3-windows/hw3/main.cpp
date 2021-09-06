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
#include "bones_loaders.hpp"
#include "bones_scene_common.hpp"


void saveScreenshot(std::string fname, void* pixels, bns::I32 w, bns::I32 h)
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

		bns::I32 _index = 1;
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
	FreeImage_Initialise(0);

	bns::I32 _index = 1;

	bns::RaytracerScene renderer_scene;
	//bns::LoadRaytracerSceneFromBnsFileFormat("triangle-ambient.bns", renderer_scene);
	bns::LoadRaytracerSceneFromBnsFileFormat("_test.test", renderer_scene);
	//bns::LoadSceneFromObjFileFormat("teapot.obj", renderer_scene);

	bns::Camera cam = renderer_scene.Cameras[0];

	bns::ColorF** colors = bns::AllocateColors(cam);
	bns::RayTrace(cam,
		renderer_scene.Shapes.data(), renderer_scene.Shapes.size(),
		renderer_scene.Lights.data(), renderer_scene.Lights.size(), colors);
	void* pixels = bns::ColorsToABGR8888Pixels(cam, colors);

	//saveScreenshot("scene_7.png", pixels, cam.ScreenWidth, cam.ScreenHeight);

	FreeImage_DeInitialise();

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