#include "bones_math.hpp"

F32 bns::IntersectionDistanceRayTriangle(const bns::RayF& ray, const bns::TriangleF& triangle)
{
	// approach: Ray-Plane intersection, then check if inside a triangle
	// ref: https://learning.edx.org/course/course-v1:UCSanDiegoX+CSE167x+2T2018/block-v1:UCSanDiegoX+CSE167x+2T2018+type@sequential+block@L9/block-v1:UCSanDiegoX+CSE167x+2T2018+type@vertical+block@vertical_9380d3229d4a
	// ref: https://www.youtube.com/watch?v=EZXz-uPyCyA

	// normal can be for any point, but let's use
	// n = ((C-A)x(B-A)) / ||(C-A)x(B-A)||
	bns::Vec3F normal = bns::Vec3F::Cross(triangle.C - triangle.A, triangle.B - triangle.A);
	normal.Normalize();

	// Plane equation, where . = dot product
	// plane = P . n - A . n = 0
	// Combine with ray equation, P_0 is ray position/origin, P_1_t is future ray position, direction with time t, so P1 is direction
	// ray = P = P_0+ P_1_t
	// (P_0 + P_1_t) . n = A . n
	// So find t
	// t = ((A . n) - (P_0 . n)) / (P_1 . n)

	// first get (P_1 . n), if 0 there is no intersection,  since ray direction is orthogonal to the normal of direction
	// which means ray direction line in the direction of a plane.
	F32 direction_dot_n = ray.Direction.Dot(normal);

	if (Abs(direction_dot_n) < 0.00001)
	{
		return MAX_F32;
	}

	// so first part ((A . n) - (P_0 . n))
	F32 t = triangle.A.Dot(normal) - (ray.Origin.Dot(normal));

	t /= direction_dot_n;

	// Find if inside a triangle parametrically ( barycentric coordinates ). Useful for texture mapping as well
	// Find alpha, beta , gamma where alpha is distance between intersection point and A
	// beta is distance between intersection point and B, and gamma between intersection point and C
	Vec3F intersection_point = ray.Origin + ray.Direction * t;

	Vec3F bc = triangle.BarycentricCoordinates(intersection_point);
	// for barycentric coordinates point is on triangle if x+y+z == 1 , but account for floating point inprecision here
	if (bc.X >= 0 && bc.X <= 1.0f
		&& bc.Y >= 0 && bc.Y <= 1.0f
		&& bc.Z >= 0 && bc.Z <= 1.0f)
	{
		if (Abs((bc.X + bc.Y + bc.Z) - 1.0f) < EPSILON)
		{
			return t;
		}
	}

	return MAX_F32;
}

void bns::IntersectionDistanceRaySphere(const RayF& ray, const SphereF& sphere, F32* out_t1, F32* out_t2)
{
	// ray = P0 + P1
	// sphere = (P-C)(P-C) - r^2 
	// so instead of P write ray
	// some algebraic manipulation should give 
	// t^2(P1.P1) + 2tP1.(P0-C) + (P0-C).(P0-C) - r^2 = 0;
	// this is quadtratic equation, so
	// at^2  + bt + c

	// ray to sphere (P0-C)
	bns::Vec3F sphere_to_ray = ray.Origin - sphere.Position;

	F32 a = ray.Direction.Dot(ray.Direction);
	F32 b = 2 * ray.Direction.Dot(sphere_to_ray);
	F32 c = sphere_to_ray.Dot(sphere_to_ray) - Pow(sphere.Radius);

	F32 discriminant = Pow(b) - 4 * a * c;

	// discriminant is less then 0, no solution
	if (discriminant < 0)
	{
		// no solution, ray has missed the sphere completly.
		*out_t1 = -1.0f;
		*out_t2 = -1.0f;
		return;
	}

	*out_t1 = (-b - Sqrt(discriminant)) / (2.0f * a);
	*out_t2 = (-b + Sqrt(discriminant)) / (2.0f * a);
}

F32 bns::RayF::IntersectionDistanceWithTriangle(const TriangleF& triangle) const
{
	return bns::IntersectionDistanceRayTriangle(*this, triangle);
}

void bns::RayF::IntersectionDistanceWithSphere(const SphereF& sphere, F32* out_t1, F32* out_t2) const 
{
	bns::IntersectionDistanceRaySphere(*this, sphere, out_t1, out_t2);
}

bns::ColorF::ColorF()
	:ColorF(0.0f, 0.0f, 0.0f)
{
}

bns::ColorF::ColorF(F32 r, F32 g, F32 b)
	: R(r), G(g), B(b), A(1.0f)
{
}

I32 bns::ColorF::ToABGR8888()
{
	U8 r = static_cast<U8>(R * 255);
	U8 g = static_cast<U8>(G * 255);
	U8 b = static_cast<U8>(B * 255);
	U8 a = static_cast<U8>(A * 255);

	I32 result = r;
	result |= g << 8;
	result |= b << 16;
	result |= a << 24;
	return result;
}

bns::Vec3F bns::TriangleF::BarycentricCoordinates(F32 x, F32 y, F32 z) const
{
	// https://www.youtube.com/watch?v=EZXz-uPyCyA
			// https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
	Vec3F v0 = B - A;
	Vec3F v1 = C - A;
	Vec3F v2 = Vec3F(x, y, z) - A;

	F32 d00 = v0.Dot(v0);
	F32 d01 = v0.Dot(v1);
	F32 d11 = v1.Dot(v1);
	F32 d20 = v2.Dot(v0);
	F32 d21 = v2.Dot(v1);

	F32 denom = d00 * d11 - d01 * d01;

	F32 b_y = (d11 * d20 - d01 * d21) / denom;
	F32 b_z = (d00 * d21 - d01 * d20) / denom;
	F32 b_x = 1.0f - b_y - b_z;

	return Vec3F(b_x, b_y, b_z);
}

bns::Vec3F bns::TriangleF::BarycentricCoordinates(const Vec3F& point) const
{
	return BarycentricCoordinates(point.X, point.Y, point.Z);
}

bns::SphereF::SphereF():
	bns::SphereF(.0f, .0f, .0f, 1.0f)
{
}
