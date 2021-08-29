#include "bones_math.hpp"
#include <iostream>

bns::Point4F bns::operator+(const Point4F& a, const Point4F& b)
{
	return Point4F(
		a.X + b.X,
		a.Y + b.Y,
		a.Z + b.Z,
		a.W + b.W
	);
}

bns::Point4F bns::operator+(const Point4F& a, const Vec4F& b)
{
	return Point4F(
		a.X + b.X,
		a.Y + b.Y,
		a.Z + b.Z,
		a.W + b.W
	);
}

bns::Point4F bns::operator*(const Point4F& a, const Vec4F& b)
{
	return Point4F(
		a.X * b.X,
		a.Y * b.Y,
		a.Z * b.Z,
		a.W * b.W
	);
}

bns::Vec4F bns::operator-(const bns::Vec4F& a, const bns::Vec4F& b)
{
	return bns::Vec4F(
		a.X - b.X,
		a.Y - b.Y,
		a.Z - b.Z,
		a.W - b.W
	);
}

bns::Vec4F bns::operator*(const bns::Vec4F& a, F32 scalar)
{
	return bns::Vec4F(
		a.X * scalar,
		a.Y * scalar,
		a.Z * scalar,
		a.W * scalar
	);
}

bns::Vec3F bns::operator+(Vec3F lhs, const Vec3F& rhs)
{
	lhs += rhs;
	return lhs;
}

bns::Vec3F bns::operator-(Vec3F lhs, const Vec3F& rhs)
{

	lhs -= rhs;
	return lhs;
}

bns::Vec3F bns::operator*(Vec3F lhs, const F32 scalar)
{
		lhs *= scalar;
		return lhs;
}

bns::Vec3F bns::operator*(F32 s, Vec3F v)
{
	v.X *= s;
	v.Y *= s;
	v.Z *= s;
	return v;
}

bns::Mat3x3F bns::operator*(F32 scalar, const bns::Mat3x3F& m)
{
	bns::Mat3x3F result(
		m.R0C0 * scalar, m.R0C1 * scalar, m.R0C2 * scalar,
		m.R1C0 * scalar, m.R1C1 * scalar, m.R1C2 * scalar,
		m.R2C0 * scalar, m.R2C1 * scalar, m.R2C2 * scalar
	);

	return result;
}

bns::Mat3x3F bns::operator+(const bns::Mat3x3F& a, const bns::Mat3x3F& b)
{
	return bns::Mat3x3F(
		a.R0C0 + b.R0C0, a.R0C1 + b.R0C1, a.R0C2 + b.R0C2,
		a.R1C0 + b.R1C0, a.R1C1 + b.R1C1, a.R1C2 + b.R1C2,
		a.R2C0 + b.R2C0, a.R2C1 + b.R2C1, a.R2C2 + b.R2C2
	);
}

bns::Vec4F bns::operator*(const Mat4x4F& m, const bns::Vec4F& v)
{
	bns::Vec4F res;

	res.X = m.R0C0 * v.X + m.R0C1 * v.Y + m.R0C2 * v.Z + m.R0C3 * v.W;

	res.Y = m.R1C0 * v.X + m.R1C1 * v.Y + m.R1C2 * v.Z + m.R1C3 * v.W;

	res.Z = m.R2C0 * v.X + m.R2C1 * v.Y + m.R2C2 * v.Z + m.R2C3 * v.W;

	res.W = m.R3C0 * v.X + m.R3C1 * v.Y + m.R3C2 * v.Z + m.R3C3 * v.W;

	return res;
}

bns::Point4F bns::operator*(const Mat4x4F& m, const Point4F& p)
{
	bns::Point4F res;

	res.X = m.R0C0 * p.X + m.R0C1 * p.Y + m.R0C2 * p.Z + m.R0C3 * p.W;

	res.Y = m.R1C0 * p.X + m.R1C1 * p.Y + m.R1C2 * p.Z + m.R1C3 * p.W;

	res.Z = m.R2C0 * p.X + m.R2C1 * p.Y + m.R2C2 * p.Z + m.R2C3 * p.W;

	res.W = m.R3C0 * p.X + m.R3C1 * p.Y + m.R3C2 * p.Z + m.R3C3 * p.W;

	return res;
}

bns::Mat4x4F bns::operator*(const Mat4x4F& a, const Mat4x4F& b)
{
	F32 R0C0 = a.R0C0 * b.R0C0 + a.R0C1 * b.R1C0 + a.R0C2 * b.R2C0 + a.R0C3 * b.R3C0;
	F32 R0C1 = a.R0C0 * b.R0C1 + a.R0C1 * b.R1C1 + a.R0C2 * b.R2C1 + a.R0C3 * b.R3C1;
	F32 R0C2 = a.R0C0 * b.R0C2 + a.R0C1 * b.R1C2 + a.R0C2 * b.R2C2 + a.R0C3 * b.R3C2;
	F32 R0C3 = a.R0C0 * b.R0C3 + a.R0C1 * b.R1C3 + a.R0C2 * b.R2C3 + a.R0C3 * b.R3C3;

	F32 R1C0 = a.R1C0 * b.R0C0 + a.R1C1 * b.R1C0 + a.R1C2 * b.R2C0 + a.R1C3 * b.R3C0;
	F32 R1C1 = a.R1C0 * b.R0C1 + a.R1C1 * b.R1C1 + a.R1C2 * b.R2C1 + a.R1C3 * b.R3C1;
	F32 R1C2 = a.R1C0 * b.R0C2 + a.R1C1 * b.R1C2 + a.R1C2 * b.R2C2 + a.R1C3 * b.R3C2;
	F32 R1C3 = a.R1C0 * b.R0C3 + a.R1C1 * b.R1C3 + a.R1C2 * b.R2C3 + a.R1C3 * b.R3C3;

	F32 R2C0 = a.R2C0 * b.R0C0 + a.R2C1 * b.R1C0 + a.R2C2 * b.R2C0 + a.R2C3 * b.R3C0;
	F32 R2C1 = a.R2C0 * b.R0C1 + a.R2C1 * b.R1C1 + a.R2C2 * b.R2C1 + a.R2C3 * b.R3C1;
	F32 R2C2 = a.R2C0 * b.R0C2 + a.R2C1 * b.R1C2 + a.R2C2 * b.R2C2 + a.R2C3 * b.R3C2;
	F32 R2C3 = a.R2C0 * b.R0C3 + a.R2C1 * b.R1C3 + a.R2C2 * b.R2C3 + a.R2C3 * b.R3C3;

	F32 R3C0 = a.R3C0 * b.R0C0 + a.R3C1 * b.R1C0 + a.R3C2 * b.R2C0 + a.R3C3 * b.R3C0;
	F32 R3C1 = a.R3C0 * b.R0C1 + a.R3C1 * b.R1C1 + a.R3C2 * b.R2C1 + a.R3C3 * b.R3C1;
	F32 R3C2 = a.R3C0 * b.R0C2 + a.R3C1 * b.R1C2 + a.R3C2 * b.R2C2 + a.R3C3 * b.R3C2;
	F32 R3C3 = a.R3C0 * b.R0C3 + a.R3C1 * b.R1C3 + a.R3C2 * b.R2C3 + a.R3C3 * b.R3C3;

	Mat4x4F result(
		R0C0, R0C1, R0C2, R0C3,
		R1C0, R1C1, R1C2, R1C3,
		R2C0, R2C1, R2C2, R2C3,
		R3C0, R3C1, R3C2, R3C3
	);
	return result;
}

std::ostream& bns::operator<<(std::ostream& os, const bns::Mat2x2F& m)
{
	os << "| " << m.R0C0 << " " << m.R0C1 << " |" << std::endl;
	os << "| " << m.R1C0 << " " << m.R1C1 << " |" << std::endl;
	return os;
}

std::ostream& bns::operator<<(std::ostream& os, const bns::Mat3x3F& m)
{
	os << "| " << m.R0C0 << " " << m.R0C1 << " " << m.R0C2 << " |" << std::endl;
	os << "| " << m.R1C0 << " " << m.R1C1 << " " << m.R1C2 << " |" << std::endl;
	os << "| " << m.R2C0 << " " << m.R2C1 << " " << m.R2C2 << " |" << std::endl;
	return os;
}

std::ostream& bns::operator<<(std::ostream& os, const bns::Mat4x4F& m)
{
	os << "| " << m.R0C0 << " " << m.R0C1 << " " << m.R0C2 << " " << m.R0C3 << " |" << std::endl;
	os << "| " << m.R1C0 << " " << m.R1C1 << " " << m.R1C2 << " " << m.R1C3 << " |" << std::endl;
	os << "| " << m.R2C0 << " " << m.R2C1 << " " << m.R2C2 << " " << m.R2C3 << " |" << std::endl;
	os << "| " << m.R3C0 << " " << m.R3C1 << " " << m.R3C2 << " " << m.R3C3 << " |" << std::endl;
	return os;
}

bns::Vec3F bns::Normal(const TriangleF& tri)
{
	Vec3F ca = tri.C - tri.A; // edge 0 
	Vec3F ba = tri.B - tri.A; // edge 1 
	Vec3F result = Vec3F::Cross(ba, ca); // this is the triangle's normal 
	result.Normalize();
	return result;
}

F32 bns::IntersectionDistanceRayTriangle(const bns::RayF& ray, const bns::TriangleF& triangle)
{
	bns::Vec3F origin = ray.Origin.ToVec3F();
	bns::Vec3F direction = ray.Direction.ToVec3F();


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
	F32 direction_dot_n = direction.Dot(normal);

	if (Abs(direction_dot_n) < 0.00001)
	{
		return MAX_F32;
	}

	// so first part ((A . n) - (P_0 . n))
	F32 t = triangle.A.Dot(normal) - (origin.Dot(normal));
	t /= direction_dot_n;

	// Find if inside a triangle parametrically ( barycentric coordinates ). Useful for texture mapping as well
	// Find alpha, beta , gamma where alpha is distance between intersection point and A
	// beta is distance between intersection point and B, and gamma between intersection point and C
	Vec3F intersection_point = (origin + direction * t);

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

	// ray to sphere (P0-C) , ignore W component here, it's not relevant
	bns::Vec3F sphere_to_ray = ray.Origin.ToVec3F() - sphere.Position;

	bns::Vec3F direction = ray.Direction.ToVec3F();

	F32 a = direction.Dot(direction);
	F32 b = 2 * direction.Dot(sphere_to_ray);
	F32 c = sphere_to_ray.Dot(sphere_to_ray) - Pow(sphere.Radius);

	F32 discriminant = Pow(b) - 4 * a * c;

	// discriminant is less then 0, no solution
	if (discriminant < 0)
	{
		// no solution, ray has missed the sphere completly.
		*out_t1 = MAX_F32;
		*out_t2 = MAX_F32;
		return;
	}

	*out_t1 = (-b - Sqrt(discriminant)) / (2.0f * a);
	*out_t2 = (-b + Sqrt(discriminant)) / (2.0f * a);
}

bns::Vec3F bns::Normalize(const bns::Vec3F& v)
{
	bns::Vec3F result = v;
	result.Normalize();
	return result;
}

F32 bns::Dot(const Vec3F& a, const Vec3F& b)
{
	F32 result = a.X * b.X + a.Y * b.Y + a.Z * b.Z;
	return result;
}

bns::RayF bns::RayF::operator*=(const bns::Mat4x4F& m)
{
	// important, notice to ToVec3, ToVec4. Point has for W by default 1, while vec has 3 when casting.
	this->Direction = m * this->Direction;
	this->Origin = m * this->Origin;

	return *this;
}

bns::RayF bns::operator*(const bns::RayF& ray, const bns::Mat4x4F& m)
{
	// important, notice to ToVec3, ToVec4. Point has for W by default 1, while vec has 3 when casting.
	bns::Vec4F dir = m * ray.Direction;
	bns::Point4F org = m * ray.Origin;

	return bns::RayF(org, dir);
}

F32 bns::RayF::IntersectionDistanceWithTriangle(const TriangleF& triangle) const
{
	return bns::IntersectionDistanceRayTriangle(*this, triangle);
}

void bns::RayF::IntersectionDistanceWithSphere(const SphereF& sphere, F32* out_t1, F32* out_t2) const
{
	bns::IntersectionDistanceRaySphere(*this, sphere, out_t1, out_t2);
}

bns::RayF::RayF(Point3F origin, Vec3F direction)
	:Origin(origin.X, origin.Y, origin.Z, 1.0f), Direction(direction.X, direction.Y, direction.Z, 0.0f)
{
}

bns::RayF::RayF(Point4F origin, Vec4F direction)
	: Origin(origin), Direction(direction)
{
}

bns::RayF::RayF(Point4F origin, Vec3F direction)
	: Origin(origin), Direction(direction.ToVec4F())
{
}

bns::ColorF::ColorF()
	: ColorF(0.0f, 0.0f, 0.0f)
{
}

bns::ColorF::ColorF(F32 r, F32 g, F32 b)
	: R(r), G(g), B(b), A(1.0f)
{
}

bns::ColorF::ColorF(F32 r, F32 g, F32 b, F32 a)
	: R(r), G(g), B(b), A(a)
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

I32 bns::ColorF::ToARGB8888()
{
	U8 r = static_cast<U8>(R * 255);
	U8 g = static_cast<U8>(G * 255);
	U8 b = static_cast<U8>(B * 255);
	U8 a = static_cast<U8>(A * 255);

	I32 result = b;
	result |= g << 8;
	result |= r << 16;
	result |= a << 24;
	return result;
}

bns::ColorF& bns::ColorF::operator+=(const ColorF& other)
{
	this->A += other.A;
	this->R += other.R;
	this->G += other.G;
	this->B += other.B;
	return *this;
}

bns::ColorF bns::ColorF::Black()
{
	return ColorF(0.0f,0.0f,0.0f,0.0f);
}

bns::ColorF bns::operator+(const bns::ColorF& a, const bns::ColorF& b)
{
	bns::ColorF result;

	result.R = a.R + b.R;
	result.G = a.G + b.G;
	result.B = a.B + b.B;
	result.A = a.A + b.A;

	return result;
}

bns::ColorF bns::operator*(const bns::ColorF& a, const bns::ColorF& b)
{
	bns::ColorF result;

	result.R = a.R * b.R;
	result.G = a.G * b.G;
	result.B = a.B * b.B;
	result.A = a.A * b.A;

	return result;
}

bns::ColorF bns::operator*(const ColorF& col, F32 scalar)
{
	return ColorF(
		col.R * scalar,
		col.G * scalar,
		col.B * scalar,
		col.A * scalar
	);
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

bns::Vec3F bns::TriangleF::GetNormal() const
{
	return Normal(*this);
}

bns::SphereF::SphereF() :
	bns::SphereF(.0f, .0f, .0f, 1.0f)
{
}

bns::Mat2x2F::Mat2x2F(
	F32 r0c0, F32 r0c1,
	F32 r1c0, F32 r1c1)
	: R0C0(r0c0), R0C1(r0c1),
	R1C0(r1c0), R1C1(r1c1)
{
}

bns::Mat2x2F::Mat2x2F()
	: Mat2x2F(
		1.0f, 0.0f,
		0.0f, 1.0f
	) {}

F32 bns::Mat2x2F::Determinant() const
{
	F32 result = this->R0C0 * this->R1C1 - this->R1C0 * this->R0C1;
	return result;
}

F32* bns::Mat2x2F::operator[](U32 index)
{
	// Get pointer to x and increase by index, then dereference
	F32* result = (&this->R0C0 + index);

	return result;
}

bns::Mat4x4F::Mat4x4F()
	: Mat4x4F(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	)
{
}

bns::Mat4x4F::Mat4x4F(F32 r0c0, F32 r0c1, F32 r0c2, F32 r0c3,
	F32 r1c0, F32 r1c1, F32 r1c2, F32 r1c3,
	F32 r2c0, F32 r2c1, F32 r2c2, F32 r2c3,
	F32 r3c0, F32 r3c1, F32 r3c2, F32 r3c3)
	:
	R0C0(r0c0), R0C1(r0c1), R0C2(r0c2), R0C3(r0c3),
	R1C0(r1c0), R1C1(r1c1), R1C2(r1c2), R1C3(r1c3),
	R2C0(r2c0), R2C1(r2c1), R2C2(r2c2), R2C3(r2c3),
	R3C0(r3c0), R3C1(r3c1), R3C2(r3c2), R3C3(r3c3)
{

}
F32* bns::Mat4x4F::operator[](U32 index)
{
	// Get pointer to x and increase by index, then dereference
	F32* result = (&this->R0C0 + index);

	return result;
}


bns::Mat4x4F bns::Mat4x4F::operator*=(const Mat4x4F& b)
{
	auto& a = *this;

	this->R0C0 = a.R0C0 * b.R0C0 + a.R0C1 * b.R1C0 + a.R0C2 * b.R2C0 + a.R0C3 * b.R3C0;
	this->R0C1 = a.R0C0 * b.R0C1 + a.R0C1 * b.R1C1 + a.R0C2 * b.R2C1 + a.R0C3 * b.R3C1;
	this->R0C2 = a.R0C0 * b.R0C2 + a.R0C1 * b.R1C2 + a.R0C2 * b.R2C2 + a.R0C3 * b.R3C2;
	this->R0C3 = a.R0C0 * b.R0C3 + a.R0C1 * b.R1C3 + a.R0C2 * b.R2C3 + a.R0C3 * b.R3C3;

	this->R1C0 = a.R1C0 * b.R0C0 + a.R1C1 * b.R1C0 + a.R1C2 * b.R2C0 + a.R1C3 * b.R3C0;
	this->R1C1 = a.R1C0 * b.R0C1 + a.R1C1 * b.R1C1 + a.R1C2 * b.R2C1 + a.R1C3 * b.R3C1;
	this->R1C2 = a.R1C0 * b.R0C2 + a.R1C1 * b.R1C2 + a.R1C2 * b.R2C2 + a.R1C3 * b.R3C2;
	this->R1C3 = a.R1C0 * b.R0C3 + a.R1C1 * b.R1C3 + a.R1C2 * b.R2C3 + a.R1C3 * b.R3C3;

	this->R2C0 = a.R2C0 * b.R0C0 + a.R2C1 * b.R1C0 + a.R2C2 * b.R2C0 + a.R2C3 * b.R3C0;
	this->R2C1 = a.R2C0 * b.R0C1 + a.R2C1 * b.R1C1 + a.R2C2 * b.R2C1 + a.R2C3 * b.R3C1;
	this->R2C2 = a.R2C0 * b.R0C2 + a.R2C1 * b.R1C2 + a.R2C2 * b.R2C2 + a.R2C3 * b.R3C2;
	this->R2C3 = a.R2C0 * b.R0C3 + a.R2C1 * b.R1C3 + a.R2C2 * b.R2C3 + a.R2C3 * b.R3C3;

	this->R3C0 = a.R3C0 * b.R0C0 + a.R3C1 * b.R1C0 + a.R3C2 * b.R2C0 + a.R3C3 * b.R3C0;
	this->R3C1 = a.R3C0 * b.R0C1 + a.R3C1 * b.R1C1 + a.R3C2 * b.R2C1 + a.R3C3 * b.R3C1;
	this->R3C2 = a.R3C0 * b.R0C2 + a.R3C1 * b.R1C2 + a.R3C2 * b.R2C2 + a.R3C3 * b.R3C2;
	this->R3C3 = a.R3C0 * b.R0C3 + a.R3C1 * b.R1C3 + a.R3C2 * b.R2C3 + a.R3C3 * b.R3C3;

	return *this;
}

bns::Mat3x3F bns::Mat4x4F::SubMatrix(U32 row, U32 col) const
{
	if (row > 3)
	{
		row = 3;
	}
	if (col > 3)
	{
		col = 3;
	}

	bns::Mat3x3F result;

	U32 i = 0;
	U32 j = 0;

	for (U32 c = 0; c < 4; c++)
	{
		if (c == col) continue;
		for (U32 r = 0; r < 4; r++)
		{
			if (r == row) continue;

			// get value from mat4x4 and move it to appropriate 3x3 position. Be careful, [] returns pointer to places.
			*result[j * 3 + i] = this->AtIndex(r, c);
			i++;
		}
		j++;
		i = 0;
	}

	return result;
}

bns::Mat2x2F bns::Mat3x3F::SubMatrix(U32 row, U32 col) const
{
	if (row > 2)
	{
		row = 2;
	}
	if (col > 2)
	{
		col = 2;
	}

	bns::Mat2x2F result;

	U32 i = 0;
	U32 j = 0;

	for (U32 c = 0; c < 3; c++)
	{
		if (c == col) continue;
		for (U32 r = 0; r < 3; r++)
		{
			if (r == row) continue;

			// get value from mat4x4 and move it to appropriate 3x3 position. Be careful, [] returns pointer to places.
			*result[j * 2 + i] = this->AtIndex(r, c);
			i++;
		}
		j++;
		i = 0;
	}

	return result;
}

bns::Mat3x3F bns::Mat3x3F::Identity()
{
	return bns::Mat3x3F();
}

F32 bns::Mat3x3F::Determinant() const
{
	F32 a = this->Cofactor(0, 0);
	F32 b = this->Cofactor(0, 1);
	F32 c = this->Cofactor(0, 2);

	F32 d = this->R0C0 * a + this->R0C1 * b + this->R0C2 * c;
	return d;
}

F32 bns::Mat3x3F::AtIndex(U32 row, U32 col) const
{
	F32 result = *(&this->R0C0 + (col * 3 + row));
	return result;
}

F32 bns::Mat4x4F::AtIndex(U32 row, U32 col) const
{
	F32 result = *(&this->R0C0 + (col * 4 + row));
	return result;
}

F32 bns::Mat3x3F::Minor(U32 row, U32 col) const
{
	bns::Mat2x2F sub = this->SubMatrix(row, col);
	F32 d = sub.Determinant();
	return d;
}

F32 bns::Mat3x3F::Cofactor(U32 row, U32 col) const
{
	I32 sign = 1;

	// For cofactor rule is
	// | + - + |
	// | - + - |
	// | + - + |

	// if odd
	if ((row + col) % 2 == 1)
	{
		sign = -1;
	}

	F32 minor = this->Minor(row, col);
	F32 result = minor * sign;
	return result;
}

bns::Mat4x4F bns::Mat4x4F::Identity()
{
	Mat4x4F result;
	return result;
}

F32 bns::Mat4x4F::Minor(U32 row, U32 col) const
{
	bns::Mat3x3F sub = this->SubMatrix(row, col);
	F32 d = sub.Determinant();
	return d;
}

bns::Mat4x4F::Mat4x4F(bns::Mat3x3F m)
	:Mat4x4F(
		m.R0C0, m.R0C1, m.R0C2, 0.0f,
		m.R1C0, m.R1C1, m.R1C2, 0.0f,
		m.R2C0, m.R2C1, m.R2C2, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	)
{
}

F32 bns::Mat4x4F::Cofactor(U32 row, U32 col) const
{
	F32 sign = 1.0f;

	// For cofactor rule is
	// | + - + - |
	// | - + - + |
	// | + - + - |
	// | - + - + |

	// if odd
	if ((row + col) % 2 == 1)
	{
		sign = -1.0f;
	}

	F32 minor = this->Minor(row, col);
	F32 result = minor * sign;

	return result;
}

F32 bns::Mat4x4F::Determinant() const
{
	F32 a = this->Cofactor(0, 0);
	F32 b = this->Cofactor(0, 1);
	F32 c = this->Cofactor(0, 2);
	F32 d = this->Cofactor(0, 3);

	F32 result = this->R0C0 * a + this->R0C1 * b + this->R0C2 * c + this->R0C3 * d;
	return result;
}

bns::Mat4x4F bns::Mat4x4F::Inverse(const bns::Mat4x4F& m)
{
	bns::Mat4x4F result = bns::Mat4x4F::Identity();

	// TODO: 
	F32 d = m.Determinant();
	for (U32 c = 0; c < 4; c++)
	{
		for (U32 r = 0; r < 4; r++)
		{
			F32 cofactor = m.Cofactor(r, c);

			// note the fliped for col in index
			*result[r * 4 + c] = cofactor / d;
		}
	}

	return result;
}

bns::Mat4x4F bns::Mat4x4F::Transpose(bns::Mat4x4F m)
{
	return Mat4x4F(
		m.R0C0, m.R1C0, m.R2C0, m.R3C0,
		m.R0C1, m.R1C1, m.R2C1, m.R3C1,
		m.R0C2, m.R1C2, m.R2C2, m.R3C2,
		m.R0C3, m.R1C3, m.R2C3, m.R3C3
	);
}

bns::Mat4x4F bns::Mat4x4F::Translate(F32 x, F32 y, F32 z)
{
	return bns::Mat4x4F(
		1.0f, 0.0f, 0.0f, x,
		0.0f, 1.0f, 0.0f, y,
		0.0f, 0.0f, 1.0f, z,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

bns::Mat4x4F bns::Mat4x4F::Translate(bns::Vec3F v)
{
	return bns::Mat4x4F::Translate(v.X, v.Y, v.Z);
}

bns::Mat4x4F bns::Mat4x4F::Scale(F32 x, F32 y, F32 z)
{
	return bns::Mat4x4F(
		x, 0.0f, 0.0f, 0.0f,
		0.0f, y, 0.0f, 0.0f,
		0.0f, 0.0f, z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

bns::Mat4x4F bns::Mat4x4F::Scale(bns::Vec3F v)
{
	return bns::Mat4x4F::Scale(v.X, v.Y, v.Z);
}

bns::Mat4x4F bns::Mat4x4F::LookAt(const Vec3F& eye, const Vec3F& center, const Vec3F& Up)
{
	// Steps
			// 1. Create a coordinate frame for the camera
			// 2. Define a rotation matrix
			// 3. Apply appropriate translation for camera ( eye ) location

			//      a          b x w
			// w = ---    u = -------       v = w x u
			//    ||a||     || b x w ||

			// a = eye - center
	Vec3F a = eye - center;
	Vec3F w = Vec3F::Normalize(a);

	Vec3F b = Vec3F::Normalize(Up);

	Vec3F b_cross_w = Vec3F::Cross(b, w);
	Vec3F b_cross_w_unit = Vec3F::Normalize(b_cross_w);

	Vec3F u = b_cross_w_unit;

	Vec3F v = Vec3F::Cross(w, u);

	// Rotation matrix
	//        | u_x u_y u_z |
	// Ruvw = | v_x v_y v_z |
	//        | w_x w_y w_z |

	// T  = -eye
	//     | R11 R12 R13 0 |  | 1 0 0 Tx |
	// M = | R21 R22 R23 0 |  | 0 1 0 Ty |
	//	   | R31 R32 R33 0 |  | 0 0 1 Tz |
	//	   | 0   0   0   1 |  | 0 0 0 1  |

	//			 | R3x3 R3x3T3x1 |
	// lookat =  | O1x3    1     |

	Mat4x4F rotation_matrix = Mat4x4F(
		u.X, u.Y, u.Z, 0,
		v.X, v.Y, v.Z, 0,
		w.X, w.Y, w.Z, 0,
		0, 0, 0, 1);

	Mat4x4F translation_matrix = Mat4x4F(
		1, 0, 0, -eye.X,
		0, 1, 0, -eye.Y,
		0, 0, 1, -eye.Z,
		0, 0, 0, 1
	);


	// You will change this return call
	Mat4x4F result = rotation_matrix * translation_matrix;
	return result;
}

bns::Mat4x4F bns::Mat4x4F::RotationMatrix(F32 theta_in_radians, bns::Vec3F axis)
{
	// I3 = identity 3x3 matrix
	//											  | xx xy xz |              | 0 -z y |
	// R(a, theta) = cos(theta)I3 + (1-cos(theta) | xy yy yz | + sin(theta) | z 0 -x | 
	//											  | xz yz zz |              | -y x 0 |
	//					a		  +            b               +           c

	axis.Normalize();

	// Unit vector components. Unit vector is v/|v|  where v is vector, and |v| is length of vector.
	F32 x = axis.X;
	F32 y = axis.Y;
	F32 z = axis.Z;

	bns::Mat3x3F a = Cos(theta_in_radians) * bns::Mat3x3F::Identity();

	bns::Mat3x3F b = (1.0f - Cos(theta_in_radians)) * bns::Mat3x3F(
		x * x, x * y, x * z,
		x * y, y * y, y * z,
		x * z, y * z, z * z);

	bns::Mat3x3F c = Sin(theta_in_radians) * bns::Mat3x3F(
		0.0f, -z, y,
		z, 0.0f, -x,
		-y, x, 0.0f);

	// You will change this return call
	bns::Mat3x3F sum = a + b + c;
	bns::Mat4x4F result(sum);
	return result;
}

bns::Vec3F& bns::Vec3F::operator+=(const Vec3F& rhs)
{
	X += rhs.X;
	Y += rhs.Y;
	Z += rhs.Z;
	return *this;
}

bns::Vec3F& bns::Vec3F::operator-=(const Vec3F& rhs)
{
	X -= rhs.X;
	Y -= rhs.Y;
	Z -= rhs.Z;
	return *this;
}

bns::Vec3F& bns::Vec3F::operator*=(const F32 scalar)
{
	X *= scalar;
	Y *= scalar;
	Z *= scalar;
	return *this;
}

bns::Vec3F& bns::Vec3F::operator/=(const F32 scalar)
{
	X /= scalar;
	Y /= scalar;
	Z /= scalar;
	return *this;
}

F32 bns::Vec3F::Length() const
{
	F32 result = X * X + Y * Y + Z * Z;
	result = Sqrt(result);
	return result;
}

void bns::Vec3F::Normalize()
{
	F32 l = Length();
	if (l > 0)
	{
		X /= l;
		Y /= l;
		Z /= l;
	}
	else
	{
		SetLengthToZero();
	}
}

void bns::Vec3F::SetLengthToZero()
{
	X = 0;
	Y = 0;
	Z = 0;
}

bns::Vec4F bns::Vec3F::ToVec4F(F32 w) const
{
	return bns::Vec4F(X, Y, Z, 0.0f);
}

bns::Point4F bns::Vec3F::ToPoint4F(F32 w) const
{
	return bns::Point4F(X, Y, Z, w);
}

bns::Vec3F bns::Vec3F::Cross(const Vec3F& a, const Vec3F& b)
{
	Vec3F result
	{
		a.Y * b.Z - b.Y * a.Z,
		a.Z * b.X - b.Z * a.X,
		a.X * b.Y - b.X * a.Y
	};
	return result;
}

inline bns::Vec3F bns::Vec3F::Normalize(const Vec3F& in)
{
	F32 l = in.Length();
	Vec3F result =
	{
		in.X / l,
		in.Y / l,
		in.Z / l,
	};
	return result;
}

inline F32 bns::Vec3F::Dot(const bns::Vec3F& a, const bns::Vec3F& b)
{
	F32 result = a.X * b.X + a.Y * b.Y + a.Z * b.Z;
	return result;
}

inline bns::Vec3F bns::Vec3F::Reflect(const Vec3F& v, const Vec3F& n)
{
	F32 v_dot_n = bns::Vec3F::Dot(v, n);
	bns::Vec3F result = v - 2 * v_dot_n * n;
	return result;
}

F32 bns::Vec3F::Dot(const Vec3F& other) const
{
	F32 result = this->X * other.X +
		this->Y * other.Y +
		this->Z * other.Z;

	return result;
}

bns::Vec3F bns::Vec3F::UnitZ()
{
	return { 0.0f, 0.0f, 1.0f };
}

bns::Vec4F::Vec4F(F32 x, F32 y, F32 z, F32 w)
	:X(x), Y(y), Z(z), W(w)
{
}

bns::Vec4F::Vec4F()
	: bns::Vec4F(0.0f, 0.0f, 0.0f, 0.0f)
{
}

F32 bns::Vec4F::Dot(const bns::Vec4F& v) const
{
	F32 result = this->X * v.X
		+ this->Y * v.Y
		+ this->Z * v.Z
		+ this->W + v.W;
	return result;
}

bns::Vec3F bns::Vec4F::ToVec3F() const
{
	return bns::Vec3F(X, Y, Z);
}

bns::Point4F::Point4F(F32 x, F32 y, F32 z, F32 w)
	: X(x), Y(y), Z(z), W(w)
{
}

bns::Point4F::Point4F()
	: Point4F(0.0f, 0.0f, 0.0f, 0.0f)
{
}

bns::Vec4F bns::Point4F::ToVec4F() const
{
	return Vec4F(X, Y, Z, W);
}

bns::Point3F bns::Point4F::ToPoint3F() const
{
	return bns::Point3F(X, Y, Z);
}

bns::Vec3F bns::Point4F::ToVec3F() const
{
	return Vec3F(X, Y, Z);
}

bns::Mat3x3F::Mat3x3F()
	:Mat3x3F(1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 1.0f)
{
}

bns::Mat3x3F::Mat3x3F(
	F32 r0c0, F32 r0c1, F32 r0c2,
	F32 r1c0, F32 r1c1, F32 r1c2,
	F32 r2c0, F32 r2c1, F32 r2c2)
	: R0C0(r0c0), R0C1(r0c1), R0C2(r0c2),
	R1C0(r1c0), R1C1(r1c1), R1C2(r1c2),
	R2C0(r2c0), R2C1(r2c1), R2C2(r2c2)
{
}

F32* bns::Mat3x3F::operator[](U32 index)
{
	// Get pointer to x and increase by index, then dereference
	F32* result = (&this->R0C0 + index);

	return result;
}

F32* bns::Vec2F::operator[](U32 index)
{
	// Get pointer to x and increase by index, then dereference
	F32* result = (&this->X + index);

	return result;
}

bns::Point3F::Point3F(F32 x, F32 y, F32 z)
	:X(x), Y(y), Z(z)
{
}

bns::Point3F::Point3F()
	: bns::Point3F(0.0f, 0.0f, 0.0f)
{
}

bns::Vec3F bns::Point3F::ToVec3F() const
{
	return bns::Vec3F(X, Y, Z);
}
