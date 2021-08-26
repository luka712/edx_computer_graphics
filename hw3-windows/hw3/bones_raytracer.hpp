#pragma once

#ifndef BONES_RAYTRACER_H

#define BONES_RAYTRACER_H

#include "bones_math.hpp"
#include "bones_material.hpp"

namespace bns
{
	/// <summary>
	/// The camera struct.
	/// </summary>
	struct Camera
	{
		bns::Vec3F LookFrom;
		bns::Vec3F LookAt;
		bns::Vec3F Up;
		U32 ScreenWidth;
		U32 ScreenHeight;
		// Field of view 
		F32 FOV;

		Camera(bns::Vec3F look_from, bns::Vec3F look_at, bns::Vec3F up, F32 fov, U32 screen_width, U32 screen_height);

		/// <summary>
		/// Get the aspect ratio.
		/// </summary>
		F32 AspectRatio() const;
		F32 FOVInRadians() const;
	};

	/// <summary>
	/// Which type of shape. Use for casting to correct type.
	/// </summary>
	enum class ShapeType
	{
		Sphere,
		Triangle
	};

	/// <summary>
	/// Just a pair of F32 values used when calcing intersection distance.
	/// </summary>
	struct F32Pair
	{
		F32 t1;
		F32 t2;
	};

	struct BaseShape
	{
		bns::Material Material;

		BaseShape();

		virtual F32Pair IntersectionDistance(const bns::RayF& ray) const = 0;
	};

	struct TriangleShape  : BaseShape
	{
		bns::TriangleF Triangle;

		TriangleShape();

		F32Pair IntersectionDistance(const bns::RayF& ray) const override;
	};

	struct SphereShape :BaseShape
	{
		bns::SphereF Sphere;

		SphereShape();

		F32Pair IntersectionDistance(const bns::RayF& ray) const override;
	};

	struct Intersection
	{
		F32 MinDist;
		const BaseShape* HitShape;
		ShapeType Type;
	};


	/// <summary>
	/// Gets a ray that passes through a pixel. 
	/// In other words gets a ray from camera to pixel on screen.
	/// </summary>
	bns::RayF RayThroughPixel(Camera cam, I32 screen_pixel_x, I32 screen_pixel_y);

	/// <summary>
	/// Find color at intersection.
	/// </summary>
	bns::ColorF FindColor(Intersection hit);

	/// <summary>
	/// Get the closest intersection between ray and triangles.
	/// 
	/// NOTE: shapes is pointer to array of BaseShape pointers. Read as array<BaseShape*>
	/// </summary>
	/// <returns></returns>
	bns::Intersection GetIntersection(const bns::RayF& ray, bns::BaseShape** shapes, U32 shapes_count);

	/// <summary>
	/// Allocate colors for raytracer.
	/// </summary>
	bns::ColorF** AllocateColors(const bns::Camera& camera);

	/// <summary>
	/// Trace the colors for shapes.
	/// 
	/// NOTE: shapes is pointer to array of BaseShape pointers. Read as array<BaseShape*>
	/// </summary>
	void RayTrace(const Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::ColorF** colors_to_fill);

	/// <summary>
	/// Get the pixels.
	/// </summary>
	void* ColorsToABGR8888Pixels(const Camera& camera, bns::ColorF** colors);

	/// <summary>
	/// Free the allocates colors.
	/// </summary>
	void FreeColors(const bns::Camera& camera, bns::ColorF** colors);

	/// <summary>
	/// Free the pixels if any were allocated.
	/// </summary>
	void FreePixels(void* pixels);
}

#endif // !BONES_RAYTRACER_H
