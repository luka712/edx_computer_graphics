#include "bones_scene_common.hpp"

#pragma region CAMERA

bns::Camera::Camera(bns::Vec3F look_from, bns::Vec3F look_at, bns::Vec3F up, bns::F32 fov_degrees, bns::U32 screen_width, bns::U32 screen_height)
	: LookFrom(look_from), LookAt(look_at), Up(up), FOVDegrees(fov_degrees), ScreenWidth(screen_width), ScreenHeight(screen_height)
{
	AspectRatio = static_cast<bns::F32>(ScreenWidth) / static_cast<bns::F32>(ScreenHeight);
	FOVRadians = bns::Radians(FOVDegrees);
}

#pragma endregion
