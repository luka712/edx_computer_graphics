#pragma once

#ifndef bones_math_H

#define bones_math_H 

#include "bones_types.hpp"
#include <math.h>
#include <cstdlib>
#include <ostream>

// NOTE: when working with intersection, if MAX_F32 is returned, there is no intersection.


// NOTE: move all of this to library, so that it can be reused between projects.

/*
	BONES MATH PROJECTS. Note to use it "bones_types.hpp" must be added.
*/

namespace bns
{
#define PI 3.14159274101257324f
#define HALF_PI (PI * 0.5f)
#define TWO_PI (PI * 2.0f)
#define EPSILON  0.0001f

	inline F32 Pow(F32 x, F32 y = 2.0f)
	{
		F32 result = powf(x, y);
		return result;
	}

	inline F32 Sqrt(F32 val)
	{
		// TODO: improve
		F32 result = sqrtf(val);
		return result;
	}

	/// <summary>
	/// Returns a random number, between two given numbers.
	/// </summary>
	/// <param name="min">Min random number.</param>
	/// <param name="max">Max random number, not including one passed in. Meaning if 4 is passed it's never going to be above 3.</param>
	/// <returns>Random number.</returns>
	inline I32 Random(I32 min, I32 max)
	{
		I32 value = max - min;
		I32 result = rand() % value + min;
		return result;
	}

	/// <summary>
	/// Returns a random number, between two given numbers.
	/// </summary>
	/// <param name="min">Min random number.</param>
	/// <param name="max">Max random number, not including one passed in. Meaning if 4 is passed it's never going to be above 3.</param>
	///// <returns>Random number.</returns>
	inline F32 Random(F32 min, F32 max)
	{
		F32 value = max - min;
		// value between 0.0f and 1.0f
		F32 random = static_cast<F32>(rand()) / static_cast<F32>(RAND_MAX);
		// move to correct space
		F32 result = random * value + min;
		return result;
	}

	/// <summary>
	/// Clamps value between min and max and returns it.
	/// </summary>
	/// <param name="min">Min value to clamp to.</param>
	/// <param name="max">Max value to clamp to.</param>
	/// <param name="val">Value to clamp.</param>
	/// <returns>Return value.</returns>
	inline F32 Clamp(F32 min, F32 max, F32 val)
	{
		if (val < min)
		{
			val = min;
			return val;
		}
		else if (val > max)
		{
			val = max;
			return val;
		}
		return val;
	}

	/// <summary>
	/// Clamps value between min and max.
	/// </summary>
	/// <param name="min">Min value to clamp to.</param>
	/// <param name="max">Max value to clamp to.</param>
	/// <param name="out_val">Value to clamp.</param>
	inline void Clamp(F32 min, F32 max, F32* out_val)
	{
		*out_val = Clamp(min, max, *out_val);
	}

	/// <summary>
	/// Get the sine of an angle value.
	/// </summary>
	/// <param name="theta">Angle theta in radians.</param>
	inline F32 Sin(F32 theta)
	{
		F32 result = sinf(theta);
		return result;
	}

	/// <summary>
	/// Get the cos of an angle value.
	/// </summary>
	/// <param name="theta">Angle theta in radians.</param>
	inline F32 Cos(F32 theta)
	{
		F32 result = cosf(theta);
		return result;
	}

	/// <summary>
	/// Gets the tan of an angle value.
	/// </summary>
	/// <param name="theta">Angle theta in radians.</param>
	/// <returns></returns>
	inline F32 Tan(F32 theta)
	{
		return tanf(theta);
	}

	/// <summary>
	/// Converts degrees to radians
	/// </summary>
	inline F32 Radians(F32 degrees)
	{
		F32 result = (degrees * PI) / 180;
		return result;
	}

	/// <summary>
	/// Converts radians to degrees 
	/// </summary>
	inline F32 Degrees(F32 radians)
	{
		F32 result = (radians * 180) / PI;
		return result;
	}

	/// <summary>
	/// Returns the value of the arc tangent of y/x, expressed in radians.
	/// </summary>
	inline F32 Atan2(F32 y, F32 x)
	{
		F32 result = atan2f(y, x);
		return result;
	}

	/// <summary>
	/// Gets the absolute value of a number.
	/// </summary>
	inline F32 Abs(F32 value)
	{
		F32 result = abs(value);
		return result;
	}

	/// <summary>
	/// Maps from range to other range.
	/// </summary>
	/// <param name="value">Amount to map from range to other range.</param>
	/// <param name="start1">Range to map from start.</param>
	/// <param name="end1">Range to map from end.</param>
	/// <param name="start2">Range to map to start.</param>
	/// <param name="end2">Range to map to end.</param>
	/// <returns>New value mapped to other range.</returns>
	inline F32 Map(F32 value, F32 start1, F32 end1, F32 start2, F32 end2)
	{
		F32 result = start2 + ((end2 - start2) / (end1 - start1)) * (value - start1);
		return result;
	}

	/// <summary>
	/// Get the max value of two.
	/// </summary>
	F32 Max(F32 a, F32 b);

	/// <summary>
	/// Get the min value of two.
	/// </summary>
	F32 Min(F32 a, F32 b);

	/// <summary>
	/// Swap the values.
	/// </summary>
	void Swap(F32* a, F32* b);
}
/*
	NOTES: magnitude/length of vector is same thing, but depends on context. If vector is line segment one can ask for it's length, if vector is representing a physical quantity
		   such as acceleration or velocity one ask for it's magnitude.
*/

namespace bns
{

	struct Vec3F;
	struct Vec4F;


	/// <summary>
	/// The structure for representing a simple point with F32 components.
	/// </summary>
	struct Point3F
	{
		F32 X;
		F32 Y;
		F32 Z;

		Point3F(F32 x, F32 y, F32 z);
		Point3F();

		/// <summary>
		/// Point to Vec3F.
		/// </summary>
		Vec3F ToVec3F() const;
	};

	/// <summary>
	/// The structure for representing a simple point with F32 components.
	/// </summary>
	struct Point4F
	{
		F32 X;
		F32 Y;
		F32 Z;
		F32 W;

		Point4F(F32 x, F32 y, F32 z, F32 w);
		Point4F();

		/// <summary>
		/// Cast point to vec 3.
		/// </summary>
		Vec3F ToVec3F() const;

		/// <summary>
		/// Cast point to vec 4.
		/// </summary>
		/// <returns></returns>
		Vec4F ToVec4F() const;

		/// <summary>
		/// Cast point to point 3F.
		/// </summary>
		Point3F ToPoint3F() const;

		friend Point4F operator+(const Point4F& a, const Point4F& b);
		friend Point4F operator+(const Point4F& a, const Vec4F& b);
		friend Point4F operator*(const Point4F& a, const Vec4F& b);
	};

	/// <summary>
	/// The Vec2 struct, which has components X and Y as F32 
	/// </summary>
	struct Vec2F
	{
		F32 X;
		F32 Y;

		/// <summary>
		/// Default constructor.
		/// Initializes properties to 0.
		/// </summary>
		Vec2F()
			: X(0.0f), Y(0.0f) {}

		/// <summary>
		/// Vec2 constructor.
		/// </summary>
		/// <param name="x">The X component.</param>
		/// <param name="y">The Y component.</param>
		Vec2F(F32 x, F32 y) :
			X(x), Y(y) {}

		/// <summary>
		/// Get the component by index.
		/// </summary>
		F32* operator[](U32 index);

		/// <summary>
		/// Set length of vector to 0.
		/// </summary>
		inline void SetLengthToZero()
		{
			X = 0;
			Y = 0;
		}

		/// <summary>
		/// Set the magnitude of vector to 0.
		/// </summary>
		inline void SetMagnitudeToZero()
		{
			SetLengthToZero();
		}

		/// <summary>
		/// Set the length of a vector.
		/// </summary>
		inline void SetLength(F32 length)
		{
			Normalize();
			X *= length;
			Y *= length;
		}

		/// <summary>
		/// Set the magnitude of a vector.
		/// </summary>
		inline void SetMagnitude(F32 magnitude)
		{
			SetLength(magnitude);
		}

		/// <summary>
		/// Get length or magnitude of a vector.
		/// </summary>
		inline F32 Length()
		{
			F32 result = X * X + Y * Y;
			result = Sqrt(result);
			return result;
		}

		/// <summary>
		/// Get the magnitude of a vector.
		/// </summary>
		/// <returns></returns>
		inline F32 Magnitude()
		{
			F32 result = Length();
			return result;
		}

		/// <summary>
		/// Gets the squared magnitude of a vector.
		/// </summary>
		/// <returns></returns>
		inline F32 MagnitudeSq()
		{
			F32 result = LengthSquared();
			return result;
		}

		/// <summary>
		/// Clamps the upper bound of length of a vector to a given length.
		/// </summary>
		/// <param name="max_length">Length to clamp to, if current length is greater then given.</param>
		inline void ClampToMaxLength(F32 max_length)
		{
			// Sqaure root is slower, therefore just use squared one.
			F32 current_length = LengthSquared();
			if (current_length > max_length * max_length)
			{
				SetLength(max_length);
			}
		}

		/// <summary>
		/// Clamps the upper bound of magniture of a vector to a given magniture.
		/// </summary>
		/// <param name="max_length">Magnitude to clamp to, if current magniture is greater then given.</param>
		inline void ClampToMaxMagnitude(F32 max_magnitude)
		{
			ClampToMaxLength(max_magnitude);
		}

		/// <summary>
		/// Normalite a vector. Sets it's length or magnitude to 1.
		/// </summary>
		inline void Normalize()
		{
			F32 l = Length();
			if (l > 0)
			{
				X /= l;
				Y /= l;
			}
			else
			{
				SetLengthToZero();
			}
		}

		/// <summary>
		/// Gets the squared length of vector.
		/// </summary>
		/// <returns>Squared length of vector.</returns>
		inline F32 LengthSquared() const
		{
			F32 result = X * X + Y * Y;
			return result;
		}

		/// <summary>
		/// Returns the angle of a vector, in radians.
		/// </summary>
		inline F32 Angle() const
		{
			F32 result = Atan2(Y, X);
			return result;
		}

		/// <summary>
		/// The direction of a vector, in radians.
		/// </summary>
		inline F32 Direction() const
		{
			F32 result = Angle();
			return result;
		}

		/// <summary>
		/// Rotate the current vector.
		/// </summary>
		/// <param name="angle">Angle to rotate by.</param>
		inline void Rotate(F32 angle)
		{
			F32 new_x = X * Cos(angle) - Y * Sin(angle);
			F32 new_y = X * Sin(angle) + Y * Cos(angle);
			X = new_x;
			Y = new_y;
		}

		/// <summary>
		/// Adds amount to current angle.
		/// </summary>
		/// <param name="angle_amount_to_add">Amount to add in radians.</param>
		inline void AddToAngle(F32 angle_amount_to_add)
		{
			F32 current_angle = Angle() + angle_amount_to_add;
			Rotate(current_angle);
		}

		/// <summary>
		/// Add amount to current direction vector.
		/// </summary>
		/// <param name="direction_amount_to_add">Amound to add in radians.</param>
		inline void AddToDirection(F32 direction_amount_to_add)
		{
			AddToAngle(direction_amount_to_add);
		}

		/// <summary>
		/// Squared distance from other vector.
		/// </summary>
		inline F32 DistanceSq(const Vec2F other)
		{
			F32 dx = this->X - other.X;
			F32 dy = this->Y - other.Y;
			F32 result = dx * dx + dy * dy;
			return result;
		}

		/// <summary>
		/// Distance from other vector.
		/// </summary>
		inline F32 Distace(const Vec2F other)
		{
			F32 result = Sqrt(DistanceSq(other));
			return result;
		}

		/// <summary>
		/// Adds vectors component wise.
		/// </summary>
		inline Vec2F& operator+=(const Vec2F& rhs)
		{
			X += rhs.X;
			Y += rhs.Y;
			return *this;
		}


		/// <summary>
		/// Adds vectors components wise.
		/// </summary>
		inline friend Vec2F operator +(Vec2F lhs, const Vec2F& rhs)
		{
			lhs += rhs;
			return lhs;
		}

		/// <summary>
		/// Subtracts vectors component wise.
		/// </summary>
		inline Vec2F& operator-=(const Vec2F& rhs)
		{
			X -= rhs.X;
			Y -= rhs.Y;
			return *this;
		}

		/// <summary>
		/// Subtracts vectors components wise.
		/// </summary>
		inline friend Vec2F operator -(Vec2F lhs, const Vec2F& rhs)
		{
			lhs -= rhs;
			return lhs;
		}

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline Vec2F& operator *=(const F32 scalar)
		{
			X *= scalar;
			Y *= scalar;
			return *this;
		}

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline friend Vec2F operator* (Vec2F lhs, const F32 scalar)
		{
			lhs *= scalar;
			return lhs;
		}

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline friend Vec2F operator* (F32 s, Vec2F v)
		{
			v.X *= s;
			v.Y *= s;
			return v;
		}

		/// <summary>
		/// Gets the zero vector.
		/// </summary>
		inline static Vec2F Zero()
		{
			Vec2F result = Vec2F(0, 0);
			return result;
		}

		/// <summary>
		/// Create the vector from angle in radians.
		/// </summary>
		inline static Vec2F FromAngle(F32 angle)
		{
			Vec2F result = { Cos(angle), Sin(angle) };
			return result;
		}

		/// <summary>
		/// Creates vector from direction in radians.
		/// </summary>
		inline static Vec2F FromDirection(F32 direction)
		{
			Vec2F result = FromAngle(direction);
			return result;
		}
	};

	/// <summary>
	/// The Vec4 struct, which has components X,Y,Z and W as F32
	/// </summary>
	struct Vec4F
	{
		F32 X;
		F32 Y;
		F32 Z;
		F32 W;

		Vec4F(F32 x, F32 y, F32 z, F32 w);
		Vec4F();

		/// <summary>
		/// Move vec4 to vec3. W is ignored.
		/// </summary>
		Vec3F ToVec3F() const ;

		/// <summary>
		/// The dot product of 2 vectors.
		/// </summary>
		F32 Dot(const bns::Vec4F& v) const;

		friend bns::Vec4F operator-(const bns::Vec4F& a, const bns::Vec4F& b);

		friend bns::Vec4F operator*(const bns::Vec4F& a, F32 scalar);
	};

	/// <summary>
	/// The Vec3 struct, which has components X,Y and Z as F32 
	/// </summary>
	struct Vec3F
	{
		F32 X;
		F32 Y;
		F32 Z;

		/// <summary>
		/// Default constructor.
		/// Initializes properties to 0.
		/// </summary>
		Vec3F()
			: X(0.0f), Y(0.0f), Z(0.0f) {}

		/// <summary>
		/// Vec2 constructor.
		/// </summary>
		/// <param name="x">The X component.</param>
		/// <param name="y">The Y component.</param>
		/// <param name="z">The Z component.</param>
		Vec3F(F32 x, F32 y, F32 z) :
			X(x), Y(y), Z(z) {}

		/// <summary>
		/// Adds vectors component wise.
		/// </summary>
		Vec3F& operator+=(const Vec3F& rhs);

		/// <summary>
		/// Adds vectors components wise.
		/// </summary>
		friend Vec3F operator +(Vec3F lhs, const Vec3F& rhs);

		/// <summary>
		/// Subtracts vectors component wise.
		/// </summary>
		Vec3F& operator-=(const Vec3F& rhs);

		/// <summary>
		/// Subtracts vectors components wise.
		/// </summary>
		friend Vec3F operator -(Vec3F lhs, const Vec3F& rhs);

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		Vec3F& operator *=(const F32 scalar);

		/// <summary>
		/// Division by scalar.
		/// </summary>
		Vec3F& operator /=(const F32 scalar);

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline friend Vec3F operator * (Vec3F lhs, const F32 scalar);

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline friend Vec3F operator* (F32 s, Vec3F v);

		/// <summary>
		/// Get length or magnitude of a vector.
		/// </summary>
		F32 Length() const;

		/// <summary>
		/// Normalite a vector. Sets it's length or magnitude to 1.
		/// </summary>
		void Normalize();

		/// <summary>
		/// Set length of vector to 0.
		/// </summary>
		void SetLengthToZero();

		/// <summary>
		/// Converts the vec3 to vec4.
		/// W components is set to 0.
		/// </summary>
		bns::Vec4F ToVec4F(F32 w = 0.0f) const;

		/// <summary>
		/// Converts the vec3 to point4.
		/// W component is set to 1.
		/// </summary>
		/// <returns></returns>
		bns::Point4F ToPoint4F(F32 w = 1.0f) const;

		/// <summary>
		/// The cross product of two vectors is the third vector that is perpendicular to the two original vectors
		/// </summary>
		/// <returns>The cross product of two vectors.</returns>
		static Vec3F Cross(const Vec3F& a, const Vec3F& b);

		/// <summary>
		/// Return the new normalized vector from in vector.
		/// </summary>
		/// <returns></returns>
		inline static Vec3F Normalize(const Vec3F& in);

		/// <summary>
		/// The dot product of two vectors.
		/// </summary>
		inline static F32 Dot(const bns::Vec3F& a, const bns::Vec3F& b);

		/// <summary>
		/// Reflect vector v around normal n.
		/// </summary>
		inline static Vec3F Reflect(const Vec3F& v, const Vec3F& n);

		/// <summary>
		/// Dot product of two vectors.
		/// a.x * b.x + a.y * b.y + a.z * b.z
		/// </summary>
		inline F32 Dot(const Vec3F& other) const;
	

		/// <summary>
		/// Retruns the vec3 with all zero components.
		/// </summary>
		inline static Vec3F Zero()
		{
			return { 0.0f, 0.0f, 0.0f };
		}

		/// <summary>
		/// Retruns the unit z vec3.
		/// </summary>
		inline static Vec3F UnitX()
		{
			return { 1.0f, 0.0f, 0.0f };
		}


		/// <summary>
		/// Retruns the unit y vec3.
		/// </summary>
		inline static Vec3F UnitY()
		{
			return { 0.0f, 1.0f, 0.0f };
		}


		/// <summary>
		/// Retruns the unit x vec3.
		/// </summary>
		static Vec3F UnitZ();
	};

	/// <summary>
	/// The color struct, with components R,G,B and A as F32
	/// </summary>
	struct ColorF
	{
		F32 R;
		F32 G;
		F32 B;
		F32 A;

		ColorF();
		ColorF(F32 r, F32 g, F32 b);
		ColorF(F32 r, F32 g, F32 b, F32 a);

		/// <summary>
		/// Return I32 where bytes are ABGR.
		/// </summary>
		I32 ToABGR8888();

		/// <summary>
		/// Return I32 where bytes are ARGB.
		/// </summary>
		I32 ToARGB8888();

		ColorF& operator+= (const ColorF & other);

		static ColorF Black();

		friend ColorF operator+(const ColorF& a, const ColorF& b);
		friend ColorF operator*(const ColorF& a, const ColorF& b);
		friend ColorF operator*(const ColorF& col, F32 scalar);
		
	};

	/// <summary>
	/// The 2x2 matrix of floats.
	/// </summary>
	struct Mat2x2F
	{

		// first col
		F32 R0C0;
		F32 R1C0;

		// second col
		F32 R0C1;
		F32 R1C1;

		Mat2x2F();

		Mat2x2F(
			F32 r0c0, F32 r0c1,
			F32 r1c0, F32 r1c1
		);

		/// <summary>
		/// Get the component by index.
		/// </summary>
		F32* operator[](U32 index);

		/// <summary>
		/// The determinant of matrix.
		/// </summary>
		F32 Determinant() const;
	};

	/// <summary>
	/// The 3x3 matrix of floats.
	/// </summary>
	struct Mat3x3F
	{
		// first col
		F32 R0C0;
		F32 R1C0;
		F32 R2C0;

		// second col
		F32 R0C1;
		F32 R1C1;
		F32 R2C1;

		// third col
		F32 R0C2;
		F32 R1C2;
		F32 R2C2;

		Mat3x3F();

		Mat3x3F(
			F32 r0c0, F32 r0c1, F32 r0c2,
			F32 r1c0, F32 r1c1, F32 r1c2,
			F32 r2c0, F32 r2c1, F32 r2c2
		);

		/// <summary>
		/// Get the component by index.
		/// </summary>
		F32* operator[](U32 index);

		/// <summary>
		/// Gets the value at index.
		/// </summary>
		F32 AtIndex(U32 row, U32 col) const;

		/// <summary>
		/// The determinant of matrix.
		/// </summary>
		F32 Determinant() const;

		/// <summary>
		/// The cofactor of a matrix.
		/// </summary>
		F32 Cofactor(U32 row, U32 col) const;

		/// <summary>
		/// The minor of a matrix.
		/// </summary>
		F32 Minor(U32 row, U32 col) const;

		/// <summary>
		/// Get the submatrix from matrix.
		/// Max value for row and col is 2 as matrices are 0 indexed.
		/// </summary>
		bns::Mat2x2F SubMatrix(U32 row, U32 col) const;

		/// <summary>
		/// Get the identity matrix.
		/// </summary>
		static bns::Mat3x3F Identity();

		/// <summary>
		/// Multiplies scalar by a matrix.
		/// </summary>
		/// <returns></returns>
		friend bns::Mat3x3F operator*(F32 scalar, const bns::Mat3x3F& m);

		/// <summary>
		/// Sums two matrices.
		/// </summary>
		friend bns::Mat3x3F operator+(const bns::Mat3x3F& a, const bns::Mat3x3F& b);
	};

#pragma region MATRIX 4x4

	/// <summary>
/// Matrix 4x4 of F32 components.
/// </summary>
	struct Mat4x4F
	{
		// first col
		F32 R0C0;
		F32 R1C0;
		F32 R2C0;
		F32 R3C0;

		// second col
		F32 R0C1;
		F32 R1C1;
		F32 R2C1;
		F32 R3C1;

		// third col
		F32 R0C2;
		F32 R1C2;
		F32 R2C2;
		F32 R3C2;

		// fourth col
		F32 R0C3;
		F32 R1C3;
		F32 R2C3;
		F32 R3C3;

		Mat4x4F();

		Mat4x4F(
			F32 r0c0, F32 r0c1, F32 r0c2, F32 r0c3,
			F32 r1c0, F32 r1c1, F32 r1c2, F32 r1c3,
			F32 r2c0, F32 r2c1, F32 r2c2, F32 r2c3,
			F32 r3c0, F32 r3c1, F32 r3c2, F32 r3c3
		);

		Mat4x4F(bns::Mat3x3F m);

		/// <summary>
		/// Get the component by index.
		/// </summary>
		F32* operator[](U32 index);

		/// <summary>
		/// Get value at matrix index.
		/// </summary>
		F32 AtIndex(U32 row, U32 col) const;

		/// <summary>
		/// Multiply matrix by other.
		/// </summary>
		Mat4x4F operator*=(const Mat4x4F& b);

		/// <summary>
		/// Multiply matrix by 3x1 vector. While this is invalid operation with matrices,
		/// here vec3 is simply casted to vec4 with w component being set to 0.
		/// </summary>
		friend Vec3F operator*(const Mat4x4F& m, const Vec3F& v);

		/// <summary>
		/// Multiply matrix by a vector.
		/// </summary>
		friend Vec4F operator*(const Mat4x4F& m, const Vec4F& v);

		/// <summary>
		/// Multiply matrix by a point.
		/// </summary>
		friend Point4F operator*(const Mat4x4F& m, const Point4F& v);


		/// <summary>
		/// Multiply two matrices.
		/// </summary>
		friend Mat4x4F operator*(const Mat4x4F& a, const Mat4x4F& b);

		/// <summary>
		/// Get the submatrix from matrix, excluding row and column.
		/// Row and col are 0 based.
		/// </summary>
		bns::Mat3x3F SubMatrix(U32 row, U32 col) const;

		/// <summary>
		/// Gets the minor of matrix.
		/// </summary>
		F32 Minor(U32 row, U32 col) const;

		/// <summary>
		/// Get the cofactor of matrix.
		/// </summary>
		F32 Cofactor(U32 row, U32 col) const;

		/// <summary>
		/// The matrix determinant.
		/// </summary>
		F32 Determinant() const;

		/// <summary>
		/// Get the identity matrix.
		/// </summary>
		static Mat4x4F Identity();

		/// <summary>
		/// The inverse of matrix.
		/// </summary>
		static Mat4x4F Inverse(const bns::Mat4x4F& m);

		/// <summary>
		/// The transpose of a matrix.
		/// </summary>
		static Mat4x4F Transpose(bns::Mat4x4F m);

		/// <summary>
		/// Creates the translation matrix.
		/// </summary>
		static Mat4x4F Translate(F32 x, F32 y, F32 z);

		/// <summary>
		/// Creates the translation matrix.
		/// </summary>
		static Mat4x4F Translate(bns::Vec3F v);

		/// <summary>
		/// Creates the scale matrix.
		/// </summary>
		static Mat4x4F Scale(F32 x, F32 y, F32 z);

		/// <summary>
		/// Creates the scale matrix.
		/// </summary>
		static Mat4x4F Scale(bns::Vec3F v);

		/// <summary>
		/// Creates rotation matrix around axis from angle theta in radians.
		/// </summary>
		static Mat4x4F RotationMatrix(F32 theta_in_radians, bns::Vec3F axis);

		/// <summary>
		/// Create model view matrix from eye ( position ), loot_at ( center ) and up vector.
		/// </summary>
		/// <returns></returns>
		static Mat4x4F LookAt(const Vec3F& eye, const Vec3F& look_at, const Vec3F& up);
	};

#pragma endregion




	/// <summary>
	/// Structure for representing a simple point with I32 components.
	/// </summary>
	struct Point2I
	{
		I32 X;
		I32 Y;
	};


	/// <summary>
	/// The rectangle with F32 components.
	/// </summary>
	struct RectF
	{
		F32 X;
		F32 Y;
		F32 W;
		F32 H;

		/// <summary>
		/// Does this rectangle intersects or overlaps other rectangle.
		/// </summary>
		/// <param name="other">The other rectangle.</param>
		/// <returns>Returns true if rectangles are intersecting.</returns>
		inline bool Intersects(const RectF& other) const
		{
			bool result = X < other.X + other.W &&
				X + W > other.X &&
				Y < other.Y + other.H &&
				Y + H > other.Y;

			return result;
		}
	};

	/// <summary>
	/// The triangle with A, B, C points
	/// </summary>
	struct TriangleF
	{
		Vec3F A;
		Vec3F B;
		Vec3F C;

		TriangleF(Vec3F a, Vec3F b, Vec3F c)
			:A(a), B(b), C(c)
		{

		}

		TriangleF()
		{
			A = bns::Vec3F::Zero();
			B = bns::Vec3F::Zero();
			C = bns::Vec3F::Zero();
		}

		/// <summary>
		/// Gets the vector that contains barycentric coordinates of a triangle from a point.
		/// If all components are <= 1 point is on triangle if x + y + z == 1
		/// </summary>
		Vec3F BarycentricCoordinates(F32 x, F32 y, F32 z) const;

		/// <summary>
		/// Gets the vector that contains barycentric coordinates of a triangle from a point.
		/// If all components are <= 1 point is on triangle if x + y + z == 1
		/// </summary>
		Vec3F BarycentricCoordinates(const Vec3F& point) const;

		/// <summary>
		/// Get the normal of triangle.
		/// </summary>
		Vec3F GetNormal() const;
	};

	Vec3F Normal(const TriangleF& tri);

	/// <summary>
	/// The sphere.
	/// </summary>
	struct SphereF
	{
		bns::Vec3F Position;
		F32 Radius;

		SphereF(F32 x, F32 y, F32 z, F32 radius)
			:Position(x, y, z), Radius(radius)
		{

		}

		SphereF();
	};


	/// <summary>
	/// The ray with origin and direction.
	/// </summary>
	struct RayF
	{
		Point4F Origin;
		Vec4F Direction;

		RayF(Point3F origin, Vec3F direction);
		RayF(Point4F origin, Vec4F direction);
		RayF(Point4F origin, Vec3F direction);

		bns::RayF operator *=(const bns::Mat4x4F& m);

		friend bns::RayF operator*(const bns::RayF& ray, const bns::Mat4x4F& m);

		/// <summary>
		/// Get the intersection distance between ray and triangle.
		/// </summary>
		/// <returns>F32, one can inspect if there is intersection by checking against MAX_F32, if result is MAX_F32 there is no intersection.</returns>
		F32 IntersectionDistanceWithTriangle(const TriangleF& triangle) const;


		/// <summary>
		/// Get the intersection distance between ray and sphere.
		/// Ray can intersect sphere at two points t1 and t2. 
		/// 
		/// t says how far does ray have to go, to hit sphere in faction P0 + P1t where P0 is ray position and P1 is ray direction:
		/// If t1 and t2 are 2 positive roots, ray is intersecting sphere. Smaller one is closer to ray.
		/// If t1 and t2 are of same positive value, ray is tangent to a sphere.
		/// If t1 and t2 are both MAX_F32 values ray has missed the sphere.
		/// If t1 or t2 is negative(complex), but other is positive value, ray is inside a sphere.
		/// </summary>
		void IntersectionDistanceWithSphere(const SphereF& sphere, F32* out_t1, F32* out_t2) const;
	};

	/// <summary>
	/// Get the intersection distance between ray and triangle.
	/// </summary>
	/// <returns>F32, one can inspect if there is intersection by checking against MAX_F32, if result is MAX_F32 there is no intersection.</returns>
	F32 IntersectionDistanceRayTriangle(const RayF& ray, const TriangleF& triangle);

	/// <summary>
	/// Get the intersection distance between ray and sphere.
	/// Ray can intersect sphere at two points t1 and t2. 
	/// 
	/// t says how far does ray have to go, to hit sphere in faction P0 + P1t where P0 is ray position and P1 is ray direction:
	/// If t1 and t2 are 2 positive roots, ray is intersecting sphere. Smaller one is closer to ray.
	/// If t1 and t2 are of same positive value, ray is tangent to a sphere.
	/// If t1 and t2 are both MAX_F32 values ray has missed the sphere.
	/// If t1 or t2 is negative(complex), but other is positive value, ray is inside a sphere.
	/// </summary>
	void IntersectionDistanceRaySphere(const RayF& ray, const SphereF& sphere, F32* out_t1, F32* out_t2);

	/// <summary>
	/// Create normalized vector, from vector being passed in.
	/// </summary>
	Vec3F Normalize(const Vec3F& v);

	/// <summary>
	/// Reflect the vector over normal.
	/// </summary>
	Vec3F Reflect(const bns::Vec3F& v, const bns::Vec3F& n);

	/// <summary>
	/// The dot product of two vectors.
	/// </summary>
	F32 Dot(const Vec3F& a, const Vec3F& b);

	/// <summary>
	/// The distance between two vectors.
	/// </summary>
	F32 Distance(const Vec3F& a, const Vec3F& b);

	std::ostream& operator<<(std::ostream& os, const Mat2x2F& m);

	std::ostream& operator<<(std::ostream& os, const Mat3x3F& m);

	std::ostream& operator<<(std::ostream& os, const Mat4x4F& m);
}
#endif 