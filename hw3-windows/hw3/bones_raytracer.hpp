#pragma once

#ifndef BONES_RAYTRACER_H

#define BONES_RAYTRACER_H

#include "bones_math.hpp"
#include "bones_material.hpp"
#include "bones_lights.hpp"

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

	struct Intersections;
	struct BaseShape;

	struct Computations
	{
		// Or min distance. It's used by ray when it intersect any of the shapes
		// P_0 + P_1T where P0 is Origin and P1 is Direction of ray.
		F32 T;
		// Intersected shape
		const BaseShape* Shape;
		// Point of intersection
		bns::Point4F LocalPoint;
		bns::Point4F WorldPoint;
		bns::Vec3F Eye;
		bns::Vec3F LocalNormal;
		bns::Vec3F WorldNormal;
	};

	struct BaseShape
	{
	private:
		bns::Mat4x4F Transform;
		bns::Mat4x4F InverseTransform;
		
	public:
		bns::Material Material;

		BaseShape();

		/// <summary>
		/// Gets the intersection distance.
		/// 
		/// NOTE: fill_intersections must be defined, and is filled by this function.
		/// </summary>
		virtual void IntersectRayAndObjects(const bns::RayF& ray, Intersections& fill_intersections) const = 0;

		/// <summary>
		/// Gets the normal
		/// </summary>
		virtual bns::Vec3F GetNormalAt(bns::Point3F point) const = 0;

		/// <summary>
		/// Gets the transform matrix.
		/// </summary>
		const bns::Mat4x4F& GetTransform() const;

		/// <summary>
		/// Gets the inverse of transform matrix. 
		/// Used like this in order to keep inverse transform cached.
		/// </summary>
		const bns::Mat4x4F& GetInverseTransform() const;

		/// <summary>
		/// Sets the transform matrix, also caches inverse transform matrix.
		/// </summary>
		void SetTransform(bns::Mat4x4F m);
	};

	struct TriangleShape  : BaseShape
	{
		bns::TriangleF Triangle;

		TriangleShape();

		void IntersectRayAndObjects(const bns::RayF& ray, Intersections& fill_intersections) const override;
		bns::Vec3F GetNormalAt(bns::Point3F point) const override;
	};

	struct SphereShape :BaseShape
	{
		bns::SphereF Sphere;

		SphereShape();

		void IntersectRayAndObjects(const bns::RayF& ray, Intersections& fill_intersections) const override;
		bns::Vec3F GetNormalAt(bns::Point3F point) const override;
	};

	struct Intersection
	{
		// TMinDist is same as t in formula P0 + P1T where P0 is Origin of ray, while P1 is Direction of ray.
		F32 TMinDist;
		const BaseShape* HitShape;
		ShapeType Type;

		Intersection operator=(const Intersection& other);
	};

	/// <summary>
	/// Struct that contains all intersections.
	/// </summary>
	struct Intersections
	{
		U32 CurrentIntersectionCount;
		// TODO: optimize this number, or maybe dynamic array ? 
		Intersection IntersectionsArray[50];

		Intersections();

		void Add(Intersection intersection);
		/// <summary>
		/// Gets the pointer intersection with minimal distance if present.
		/// </summary>
		Intersection* Hit();
	};

	/// <summary>
	/// Gets a ray that passes through a pixel. 
	/// In other words gets a ray from camera to pixel on screen.
	/// </summary>
	bns::RayF RayThroughPixel(Camera cam, I32 screen_pixel_x, I32 screen_pixel_y);

	/// <summary>
	/// Find color at intersection.
	/// </summary>
	bns::ColorF ColorAt(Intersection hit);

	/// <summary>
	/// Get the closest intersection between ray and triangles.
	/// 
	/// NOTE: shapes is pointer to array of BaseShape pointers. Read as array<BaseShape*>
	bns::Intersections GetIntersections(const bns::RayF& ray, bns::BaseShape** Shapes, U32 shapes_count);

	/// <summary>
	/// Allocate colors for raytracer.
	/// </summary>
	bns::ColorF** AllocateColors(const bns::Camera& camera);

	/// <summary>
	/// Trace the colors for shapes.
	/// 
	/// NOTE: shapes is pointer to array of BaseShape pointers. Read as array<BaseShape*>
	/// </summary>
	void RayTrace(const Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::BaseLight** lights, U32 count_of_lights, bns::ColorF** colors_to_fill);

	/// <summary>
	/// Prepare all the computations necessary for further processing of intersections, stuff as light, shadows, reflection etc...
	/// </summary>
	bns::Computations PrepareComputations(const Intersection& intersection, const RayF& ray, const Intersections& other);

	/// <summary>
	/// Applies the Blinn-Phong reflection model lighting for ray and it's computations.
	/// </summary>
	bns::ColorF BlinnPhongReflectionModel( const Computations& comp, bns::BaseShape** shapes, U32 count_of_shapes, BaseLight** lights, U32 count_of_lights);

	bns::ColorF ReflectedColor(const Computations& comp, I32 remaining = 5);

	bool IsShadowed(const Computations& comp, const BaseLight& light, bns::BaseShape** shapes, U32 count_of_shapes);

	/// <summary>
	/// Get the pixels in format ABGR8888.
	/// </summary>
	void* ColorsToABGR8888Pixels(const Camera& camera, bns::ColorF** colors);

	/// <summary>
	/// Get the pixels in format ABGR8888.
	/// </summary>
	void* ColorsToARGB8888Pixels(const Camera& camera, bns::ColorF** colors);

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
