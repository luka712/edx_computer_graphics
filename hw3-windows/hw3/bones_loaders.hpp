#pragma once

#ifndef BONES_LOADERS_H

#define BONES_LOADERS_H

#include "bones_raytracer.hpp"

namespace bns
{
	/// <summary>
	/// Load the scene from .bns file format which is native to bns library.
	/// </summary>
	void LoadRaytracerSceneFromBnsFileFormat(const char* filename, RaytracerScene& scene);

	/// <summary>
	/// Load the scene from .obj file format ( wavefront)
	/// </summary>
	void LoadSceneFromObjFileFormat(const char* filename, RaytracerScene& scene);

}

#endif // !BONES_LOADERS_H
