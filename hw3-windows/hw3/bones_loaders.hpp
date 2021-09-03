#pragma once

#ifndef BONES_LOADERS_H

#define BONES_LOADERS_H

#include "bones_raytracer.hpp"

namespace bns
{
	/// <summary>
	/// Load the scene from .bns file format which is native to bns library.
	/// </summary>
	void LoadSceneFromBnsFileFormat(const char* filename, RayTracerScene& scene);
}

#endif // !BONES_LOADERS_H
