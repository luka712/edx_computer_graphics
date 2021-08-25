#pragma once

#ifndef bones_math_H

#define bones_math_H 

#include "bones_types.hpp"
#include <math.h>
#include <cstdlib>
#include <ostream>


// NOTE: move all of this to library, so that it can be reused between projects.

/*
	BONES MATH PROJECTS. Note to use it "bones_types.hpp" must be added.
*/

namespace bns
{
#define PI 3.14159274101257324f
#define HALF_PI (PI * 0.5f)
#define TWO_PI (PI * 2.0f)

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
}
/*
	NOTES: magnitude/length of vector is same thing, but depends on context. If vector is line segment one can ask for it's length, if vector is representing a physical quantity
		   such as acceleration or velocity one ask for it's magnitude.
*/

namespace bns
{
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
		inline F32 operator[](U32 index)
		{
			// Get pointer to x and increase by index, then dereference
			F32 result = *(&this->X + index);

			return result;
		}

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
		inline Vec3F& operator+=(const Vec3F& rhs)
		{
			X += rhs.X;
			Y += rhs.Y;
			Z += rhs.Z;
			return *this;
		}


		/// <summary>
		/// Adds vectors components wise.
		/// </summary>
		inline friend Vec3F operator +(Vec3F lhs, const Vec3F& rhs)
		{
			lhs += rhs;
			return lhs;
		}

		/// <summary>
		/// Subtracts vectors component wise.
		/// </summary>
		inline Vec3F& operator-=(const Vec3F& rhs)
		{
			X -= rhs.X;
			Y -= rhs.Y;
			Z -= rhs.Z;
			return *this;
		}

		/// <summary>
		/// Subtracts vectors components wise.
		/// </summary>
		inline friend Vec3F operator -(Vec3F lhs, const Vec3F& rhs)
		{
			lhs -= rhs;
			return lhs;
		}

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline Vec3F& operator *=(const F32 scalar)
		{
			X *= scalar;
			Y *= scalar;
			Z *= scalar;
			return *this;
		}

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline friend Vec3F operator* (Vec3F lhs, const F32 scalar)
		{
			lhs *= scalar;
			return lhs;
		}

		/// <summary>
		/// Multiplication by scalar.
		/// </summary>
		inline friend Vec3F operator* (F32 s, Vec3F v)
		{
			v.X *= s;
			v.Y *= s;
			v.Z *= s;
			return v;
		}

		/// <summary>
		/// Get length or magnitude of a vector.
		/// </summary>
		inline F32 Length() const
		{
			F32 result = X * X + Y * Y + Z * Z;
			result = Sqrt(result);
			return result;
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
				Z /= l;
			}
			else
			{
				SetLengthToZero();
			}
		}

		/// <summary>
		/// Set length of vector to 0.
		/// </summary>
		inline void SetLengthToZero()
		{
			X = 0;
			Y = 0;
			Z = 0;
		}

		/// <summary>
		/// The cross product of two vectors is the third vector that is perpendicular to the two original vectors
		/// </summary>
		/// <returns>The cross product of two vectors.</returns>
		inline static Vec3F Cross(const Vec3F& a, const Vec3F& b)
		{
			Vec3F result
			{
				a.Y * b.Z - b.Y * a.Z,
				a.Z * b.X - b.Z * a.X,
				a.X * b.Y - b.X * a.Y
			};
			return result;
		}

		/// <summary>
		/// Return the new normalized vector from in vector.
		/// </summary>
		/// <returns></returns>
		inline static Vec3F Normalize(const Vec3F& in)
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

		/// <summary>
		/// Dot product of two vectors.
		/// a.x * b.x + a.y * b.y + a.z * b.z
		/// </summary>
		inline F32 Dot(const Vec3F& other) const
		{
			F32 result = this->X * other.X +
				this->Y * other.Y +
				this->Z * other.Z;

			return result;
		}

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
		inline static Vec3F UnitZ()
		{
			return { 0.0f, 0.0f, 1.0f };
		}
	};

	struct Mat4x4F
	{
		// first column
		F32 r0c0;
		F32 r1c0;
		F32 r2c0;
		F32 r3c0;

		// second column
		F32 r0c1;
		F32 r1c1;
		F32 r2c1;
		F32 r3c1;

		// third column
		F32 r0c2;
		F32 r1c2;
		F32 r2c2;
		F32 r3c2;

		// fourth column
		F32 r0c3;
		F32 r1c3;
		F32 r2c3;
		F32 r3c3;

		Mat4x4F() :Mat4x4F(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		)
		{

		}

		Mat4x4F(
			F32 r0c0, F32 r0c1, F32 r0c2, F32 r0c3,
			F32 r1c0, F32 r1c1, F32 r1c2, F32 r1c3,
			F32 r2c0, F32 r2c1, F32 r2c2, F32 r2c3,
			F32 r3c0, F32 r3c1, F32 r3c2, F32 r3c3
		) :
			r0c0(r0c0), r0c1(r0c1), r0c2(r0c2), r0c3(r0c3),
			r1c0(r1c0), r1c1(r1c1), r1c2(r1c2), r1c3(r1c3),
			r2c0(r2c0), r2c1(r2c1), r2c2(r2c2), r2c3(r2c3),
			r3c0(r3c0), r3c1(r3c1), r3c2(r3c2), r3c3(r3c3)
		{

		}

		inline friend Mat4x4F operator*(const Mat4x4F& a, const Mat4x4F& b)
		{
			F32 r0c0 = a.r0c0 * b.r0c0 + a.r0c1 * b.r1c0 + a.r0c2 * b.r2c0 + a.r0c3 * b.r3c0;
			F32 r0c1 = a.r0c0 * b.r0c1 + a.r0c1 * b.r1c1 + a.r0c2 * b.r2c1 + a.r0c3 * b.r3c1;
			F32 r0c2 = a.r0c0 * b.r0c2 + a.r0c1 * b.r1c2 + a.r0c2 * b.r2c2 + a.r0c3 * b.r3c2;
			F32 r0c3 = a.r0c0 * b.r0c3 + a.r0c1 * b.r1c3 + a.r0c2 * b.r2c3 + a.r0c3 * b.r3c3;

			F32 r1c0 = a.r1c0 * b.r0c0 + a.r1c1 * b.r1c0 + a.r1c2 * b.r2c0 + a.r1c3 * b.r3c0;
			F32 r1c1 = a.r1c0 * b.r0c1 + a.r1c1 * b.r1c1 + a.r1c2 * b.r2c1 + a.r1c3 * b.r3c1;
			F32 r1c2 = a.r1c0 * b.r0c2 + a.r1c1 * b.r1c2 + a.r1c2 * b.r2c2 + a.r1c3 * b.r3c2;
			F32 r1c3 = a.r1c0 * b.r0c3 + a.r1c1 * b.r1c3 + a.r1c2 * b.r2c3 + a.r1c3 * b.r3c3;

			F32 r2c0 = a.r2c0 * b.r0c0 + a.r2c1 * b.r1c0 + a.r2c2 * b.r2c0 + a.r2c3 * b.r3c0;
			F32 r2c1 = a.r2c0 * b.r0c1 + a.r2c1 * b.r1c1 + a.r2c2 * b.r2c1 + a.r2c3 * b.r3c1;
			F32 r2c2 = a.r2c0 * b.r0c2 + a.r2c1 * b.r1c2 + a.r2c2 * b.r2c2 + a.r2c3 * b.r3c2;
			F32 r2c3 = a.r2c0 * b.r0c3 + a.r2c1 * b.r1c3 + a.r2c2 * b.r2c3 + a.r2c3 * b.r3c3;

			F32 r3c0 = a.r3c0 * b.r0c0 + a.r3c1 * b.r1c0 + a.r3c2 * b.r2c0 + a.r3c3 * b.r3c0;
			F32 r3c1 = a.r3c0 * b.r0c1 + a.r3c1 * b.r1c1 + a.r3c2 * b.r2c1 + a.r3c3 * b.r3c1;
			F32 r3c2 = a.r3c0 * b.r0c2 + a.r3c1 * b.r1c2 + a.r3c2 * b.r2c2 + a.r3c3 * b.r3c2;
			F32 r3c3 = a.r3c0 * b.r0c3 + a.r3c1 * b.r1c3 + a.r3c2 * b.r2c3 + a.r3c3 * b.r3c3;

			Mat4x4F result(
				r0c0, r0c1, r0c2, r0c3,
				r1c0, r1c1, r1c2, r1c3,
				r2c0, r2c1, r2c2, r2c3,
				r3c0, r3c1, r3c2, r3c3
			);
			return result;
		}

		inline static Mat4x4F Identity()
		{
			Mat4x4F result;
			return result;
		}

		inline static Mat4x4F LookAt(const Vec3F& eye, const Vec3F& center, const Vec3F& up)
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

			Vec3F b = Vec3F::Normalize(up);

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

		inline friend std::ostream& operator<<(std::ostream& os, const Mat4x4F& m)
		{
			os << "| " << m.r0c0 << " " << m.r0c1 << " " << m.r0c2 << " " << m.r0c3 << " |" << std::endl;
			os << "| " << m.r1c0 << " " << m.r1c1 << " " << m.r1c2 << " " << m.r1c3 << " |" << std::endl;
			os << "| " << m.r2c0 << " " << m.r2c1 << " " << m.r2c2 << " " << m.r2c3 << " |" << std::endl;
			os << "| " << m.r3c0 << " " << m.r3c1 << " " << m.r3c2 << " " << m.r3c3 << " |" << std::endl;
			return os;
		}
	};

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

		/// <summary>
		/// Gets the vector that contains barycentric coordinates of a triangle from a point.
		/// If all components are <= 1 point is on triangle if x + y + z == 1
		/// </summary>
		inline Vec3F BarycentricCoordinates(F32 x, F32 y, F32 z) const 
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

		Vec3F BarycentricCoordinates(const Vec3F& point) const 
		{
			return BarycentricCoordinates(point.X, point.Y, point.Z);
		}
	};


	/// <summary>
	/// The ray with origin and direction.
	/// </summary>
	struct RayF
	{
		Vec3F Origin;
		Vec3F Direction;

		RayF(Vec3F origin, Vec3F direction)
			:Origin(origin), Direction(direction) {}

		inline F32 Intersects(const TriangleF& triangle) const 
		{
			// approach: Ray-Plane intersection, then check if inside a triangle
			// ref: https://learning.edx.org/course/course-v1:UCSanDiegoX+CSE167x+2T2018/block-v1:UCSanDiegoX+CSE167x+2T2018+type@sequential+block@L9/block-v1:UCSanDiegoX+CSE167x+2T2018+type@vertical+block@vertical_9380d3229d4a
			// ref: https://www.youtube.com/watch?v=EZXz-uPyCyA

			// normal can be for any point, but let's use
			// n = ((C-A)x(B-A)) / ||(C-A)x(B-A)||
			Vec3F normal = Vec3F::Cross(triangle.C - triangle.A, triangle.B - triangle.A);
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
			F32 direction_dot_n = Direction.Dot(normal);

			if (Abs(direction_dot_n) < 0.00001)
			{
				return MAX_F32;
			}

			// so first part ((A . n) - (P_0 . n))
			F32 t = triangle.A.Dot(normal) - (Origin.Dot(normal));

			t /= direction_dot_n;

			// Find if inside a triangle parametrically ( barycentric coordinates ). Useful for texture mapping as well
			// Find alpha, beta , gamma where alpha is distance between intersection point and A
			// beta is distance between intersection point and B, and gamma between intersection point and C
			Vec3F intersection_point = Origin + Direction * t;
	
			Vec3F bc = triangle.BarycentricCoordinates(intersection_point);
			// for barycentric coordinates point is on triangle if x+y+z == 1 , but account for floating point inprecision here
			if (bc.X >= 0 && bc.X <= 1.0f 
				&& bc.Y >= 0 && bc.Y <= 1.0f 
				&& bc.Z >= 0 && bc.Z <= 1.0f) 
			{
				if (Abs((bc.X + bc.Y + bc.Z) - 1.0f) < 0.00001f)
				{
					return t;
				}
			}

			return MAX_F32;
		}
	};



}
#endif 