#pragma once

#ifndef BONES_LIGHTS_H

#define BONES_LIGHTS_H

#include "bones_math.hpp"

namespace bns 
{
	enum LightType
	{
		Directional,
		Point 
	};

	struct BaseLight
	{
		bns::ColorF Color;
		bns::LightType Type;

		BaseLight(LightType type);
	};

	struct DirectionalLight final : BaseLight
	{
		bns::Vec3F Direction;
		DirectionalLight();
	};

	struct PointLight final : BaseLight
	{
		bns::Vec3F Position;
		PointLight();
	};
}


#endif // !BONES_LIGHTS_H

