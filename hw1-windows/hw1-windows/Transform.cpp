// Transform.cpp: implementation of the Transform class.

#include "Transform.h"

//Please implement the following functions:
float length(vec3 v)
{
	float result = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	return result;
}

vec3 normal_vector(vec3 v)
{
	float l = length(v);
	vec3 result = vec3(v.x / l, v.y / l, v.z / l);
	return result;
}

// Helper rotation function.  
mat3 Transform::rotate(const float degrees, const vec3& axis) {
	// YOUR CODE FOR HW1 HERE

	// theta in radians
	const float theta = degrees * pi / 180.0f;

	// I3 = identity 3x3 matrix
	//											  | xx xy xz |              | 0 -z y |
	// R(a, theta) = cos(theta)I3 + (1-cos(theta) | xy yy yz | + sin(theta) | z 0 -x | 
	//											  | xz yz zz |              | -y x 0 |
	//					a		  +            b               +           c

	float l = length(axis);

	// Unit vector components. Unit vector is v/|v|  where v is vector, and |v| is length of vector.
	float x = axis.x / l;
	float y = axis.y / l;
	float z = axis.z / l;

	mat3 a = cos(theta) * mat3(1);
	mat3 b = (1.0f - cos(theta)) * glm::transpose(mat3(x * x, x * y, x * z, x * y, y * y, y * z, x * z, y * z, z * z));
	mat3 c = sin(theta) * glm::transpose(mat3(0, -z, y, z, 0, -x, -y, x, 0));

	// You will change this return call
	mat3 result = a + b + c;
	return result;
}

// Transforms the camera left around the "crystal ball" interface
void Transform::left(float degrees, vec3& eye, vec3& up) {
	eye = rotate(degrees, up) * eye;
}

// Transforms the camera up around the "crystal ball" interface
void Transform::up(float degrees, vec3& eye, vec3& up) {
	// YOUR CODE FOR HW1 HERE 
	vec3 right = glm::cross(eye, up);
	eye = rotate(degrees, right) * eye;
	up = glm::cross(glm::normalize(right),glm::normalize(eye));
}

// Your implementation of the glm::lookAt matrix
mat4 Transform::lookAt(vec3 eye, vec3 up) {
	// YOUR CODE FOR HW1 HERE

	// Steps
	// 1. Create a coordinate frame for the camera
	// 2. Define a rotation matrix
	// 3. Apply appropriate translation for camera ( eye ) location

	// Center in opengl
	vec3 center = vec3(0, 0, 0);

	//      a          b x w
	// w = ---    u = -------       v = w x u
	//    ||a||     || b x w ||

	// a = eye - center
	vec3 a = eye - center;
	float l = length(a);
	vec3 w = vec3(a.x / l, a.y / l, a.z / l);

	l = length(up);
	vec3 b = vec3(up.x / l, up.y / l, up.z / l);


	vec3 b_cross_w = glm::cross(b, w);
	l = length(b_cross_w);
	vec3 b_cross_w_unit = vec3(b_cross_w.x / l, b_cross_w.y / l, b_cross_w.z / l);

	vec3 u = b_cross_w_unit;

	vec3 v = glm::cross(w, u);

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

	mat4 rotation_matrix = glm::transpose(mat4(
		u.x, u.y, u.z, 0,
		v.x, v.y, v.z, 0,
		w.x, w.y, w.z, 0,
		0, 0, 0, 1));

	mat4 translation_matrix = glm::transpose(mat4(
		1, 0, 0, -eye.x,
		0, 1, 0, -eye.y,
		0, 0, 1, -eye.z,
		0, 0, 0, 1
	));


	// You will change this return call
	mat4 result = rotation_matrix * translation_matrix;
	return result;
}

Transform::Transform()
{

}

Transform::~Transform()
{

}