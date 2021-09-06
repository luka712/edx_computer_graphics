#pragma once

#ifndef BONES_SCENE_COMMON_H

#define BONES_SCENE_COMMON_H

#include "bones_math.hpp"

namespace bns
{
	/// <summary>
	/// The camera struct.
	/// </summary>
	struct Camera
	{
		Vec3F LookFrom;
		Vec3F LookAt;
		Vec3F Up;

		U32 ScreenWidth;
		U32 ScreenHeight;

		F32 AspectRatio;

		// Field of view 
		F32 FOVDegrees;
		F32 FOVRadians;

		Camera(bns::Vec3F look_from, bns::Vec3F look_at, bns::Vec3F up, F32 fov_degrees, U32 screen_width, U32 screen_height);
	};
}

#endif // !BONES_SCENE_COMMON_H
