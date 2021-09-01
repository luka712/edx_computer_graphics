#pragma once
#include "bones_math.hpp"

#ifndef BONES_MATERIAL_H

#define BONES_MATERIAL_H


namespace bns
{
	struct Material
	{
		bns::ColorF Ambient;
		bns::ColorF Emission;
		bns::ColorF Diffuse;
		bns::ColorF Specular;
		F32 Shininess;
		F32 RefractiveIndex;
	};
}

#endif // !BONES_MATERIAL_H
