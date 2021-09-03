#pragma once

#ifndef BONES_RAYTRACER_H

#define BONES_RAYTRACER_H

#include <vector>
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
		Vec3F LookFrom;
		Vec3F LookAt;
		Vec3F Up;
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
		Triangle,
		Cube,
		Plane, 
		Cylinder
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

		bns::Point4F RayOrigin;
		bns::Vec3F RayDirection;

		// Point of intersection
		bns::Point4F LocalPoint;
		bns::Point4F WorldPoint;
		
		bns::Vec3F Eye;
		bns::Vec3F LocalNormal;
		bns::Vec3F WorldNormal;
		bns::Vec3F ReflectedVector;

		// refractive index of the material that ray is passing from ( being exited from )
		F32 N1;
		// refractive index of the material that ray is passing to ( being entered to )
		F32 N2;
	};



#pragma region BASE SHAPE

	struct BaseShape
	{
	private:
		bns::Mat4x4F Transform;
		bns::Mat4x4F InverseTransform;

	public:
		bns::Material Material;
		bns::BaseShape* Parent;

		BaseShape();

		/// <summary>
		/// Gets the intersection distance.
		/// 
		/// NOTE: fill_intersections must be defined, and is filled by this function.
		/// </summary>
		virtual void IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const = 0;

		/// <summary>
		/// Gets the normal
		/// </summary>
		virtual bns::Vec3F GetLocalNormalAt(bns::Point3F point) const = 0;

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

#pragma endregion

#pragma region TRIANGLE SHAPE
	
	/// <summary>
	/// The triangle shape.
	/// </summary>
	struct TriangleShape : BaseShape
	{
		bns::TriangleF Triangle;

		TriangleShape();

		void IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const override;
		bns::Vec3F GetLocalNormalAt(bns::Point3F point) const override;
	};

#pragma endregion

#pragma region SPHERE SHAPE

	struct SphereShape :BaseShape
	{
		bns::SphereF Sphere;

		SphereShape();

		/// <summary>
		/// Checks if there are intersections between ray and shape and fills intersections accordingly.
		/// </summary>
		void IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const override;

		/// <summary>
		/// Gets the local normal of shape at passed point.
		/// </summary>
		bns::Vec3F GetLocalNormalAt(bns::Point3F point) const override;
	};

#pragma endregion

#pragma region CUBE SHAPE

	struct CubeShape : BaseShape
	{
		CubeShape();

		/// <summary>
		/// Checks if ray intersects plane and return minimum and maximum t value.
		/// </summary>
		void CheckAxis(F32 origin, F32 direction, F32* tmin, F32* tmax) const;
		void IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const override;
		bns::Vec3F GetLocalNormalAt(bns::Point3F point) const override;
	};

#pragma endregion

#pragma region PLANE SHAPE

	struct PlaneShape : BaseShape
	{
		void IntersectWithRay(const RayF& ray, Intersections& fill_intersections) const override;
		Vec3F GetLocalNormalAt(bns::Point3F point) const override;
	};

#pragma endregion

#pragma region CYLINDER SHAPE

	struct CylinderShape : BaseShape
	{
		// Minimum and Maximum define edges of cylinder in order for it, not to stretch to infinity.
		F32 Minimum;
		F32 Maximum;

		// Is cylinder capped at bottom and top.
		bool Capped;

		CylinderShape();
		void IntersectWithRay(const RayF& ray, Intersections& fill_intersections) const override;
		Vec3F GetLocalNormalAt(bns::Point3F point) const override;

	private:
		/// <summary>
		/// Helper function.
		/// Checks to see if the intersection as `t` is withing a radius
		/// of 1 ( the radius of cylinders ) from the y axis.
		/// </summary>
		bool CheckCap(const RayF& ray, F32 t) const;

		/// <summary>
		/// Helper function.
		/// Fills intersections for cylinders caps.
		/// </summary>
		/// <param name="ray"></param>
		/// <param name="fill_intersections"></param>
		void IntersectCaps(const RayF& ray, Intersections& fill_intersections) const;
	};

#pragma endregion

#pragma region GROUP SHAPE

	class Group : public BaseShape
	{

	};

#pragma endregion

#pragma region SHAPE HELPERS

	/// <summary>
	/// Transform the local_normal to world normal.
	/// </summary>
	bns::Vec3F GetWorldNormalAt(const bns::BaseShape& shape, bns::Vec3F local_normal);

	/// <summary>
	/// Transform the world point to local point.
	/// </summary>
	bns::Point4F WorldPointToLocalPoint(const bns::BaseShape& shape, bns::Point4F world_point);

#pragma endregion

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
		// TODO: Move these to heap!!!!
		Intersection IntersectionsArray[10000];

		Intersections();

		void Add(Intersection intersection);
		/// <summary>
		/// Gets the pointer intersection with minimal distance if present.
		/// </summary>
		Intersection* Hit();
	};

	/// <summary>
	/// The raytracer scene.
	/// </summary>
	struct RayTracerScene
	{
		std::vector<Camera> Cameras;
		std::vector<BaseShape*> Shapes;
		std::vector<BaseLight*> Lights;
	};

	/// <summary>
	/// Gets a ray that passes through a pixel. 
	/// In other words gets a ray from camera to pixel on screen.
	/// </summary>
	bns::RayF RayThroughPixel(Camera cam, I32 screen_pixel_x, I32 screen_pixel_y);

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
	/// Get the color from ray.
	/// </summary>
	bns::ColorF ColorAt(const bns::RayF& ray,
		bns::BaseShape** shapes, U32 count_of_shapes,
		bns::BaseLight** lights, U32 count_of_lights,
		I32 remaining);

	/// <summary>
	/// Trace the colors for shapes, uses multiple threads to work.
	/// 
	/// NOTE: shapes is pointer to array of BaseShape pointers. Read as array<BaseShape*>
	/// </summary>
	void ThreadedRayTrace(const Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::BaseLight** lights, U32 count_of_lights, bns::ColorF** colors_to_fill);

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

	/// <summary>
	/// Get the reflected color.
	/// 
	/// This function created the reflection ray and calls ColorAt for it, in other words gets it's lighting, shadow and further calls reflected rays for it.
	/// </summary>
	bns::ColorF ReflectedColor(const Computations& comp, 
		bns::BaseShape** shapes, U32 count_of_shapes,
		bns::BaseLight** lights, U32 count_of_lights,
		I32 remaining);

	/// <summary>
	/// Get the refracted color.
	/// 
	/// This function created the refraction ray and calls ColorAt for it, in other words gets it's lighting, shadow and further calls refracted rays for it.
	/// </summary>
	bns::ColorF RefractedColor(const Computations& comp,
		bns::BaseShape** shapes, U32 count_of_shapes,
		bns::BaseLight** lights, U32 count_of_lights,
		I32 remaining);

	/// <summary>
	/// Returns boolean expressions which says if point in space is directly blocked by other object from light source.
	/// Determines if object is in shadow or not.
	/// </summary>
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
