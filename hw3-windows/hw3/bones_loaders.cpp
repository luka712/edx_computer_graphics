#include "bones_loaders.hpp"
#include <vector>
#include <sstream>
#include <stack>
#include "bones_math.hpp"
#include "bones_files.hpp"
#include "bones_char_string.hpp"


void bns::LoadSceneFromBnsFileFormat(const char* filename, bns::RayTracerScene& scene)
{
	std::vector<bns::Vec3F> vertices;

	bns::I32 camera_width = 0;
	bns::I32 camera_height = 0;

	bns::FileContents* file = bns::ReadAndCloseFile(filename);

	std::string line;
	std::stringstream file_stream(file->Contents);


	std::stack<bns::Mat4x4F> matrix_stack;
	matrix_stack.push(bns::Mat4x4F::Identity());

	bns::ColorF ambient_color = bns::ColorF(.2f, .2f, .2f);
	bns::ColorF emission_color;
	bns::ColorF diffuse_color;
	bns::ColorF specular_color;
	bns::F32 shininess = 0.0f;
	bns::F32 attenuation_constant = 1.0f;
	bns::F32 attenuation_linear = 0.0f;
	bns::F32 attenuation_quadratic = 0.0f;

	while (std::getline(file_stream, line))
	{
		if (line.rfind("end", 0) == 0)
		{
			break;
		}
		if (line.rfind("size", 0) == 0)
		{
			bns::I32 end_index = 0;

			bns::ReadI32FromString(line.c_str(), end_index, &end_index, &camera_width);
			bns::ReadI32FromString(line.c_str(), end_index, &end_index, &camera_height);
		}
		else if (line.rfind("camera", 0) == 0)
		{
			bns::I32 end_index = 0;
			const char* ptr_to_camera = line.c_str();

			bns::Vec3F eyeinit;
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &eyeinit.X);
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &eyeinit.Y);
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &eyeinit.Z);

			bns::Vec3F center;
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &center.X);
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &center.Y);
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &center.Z);

			bns::Vec3F upinit;
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &upinit.X);
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &upinit.Y);
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &upinit.Z);

			bns::F32 fov;
			bns::ReadF32FromString(ptr_to_camera, end_index, &end_index, &fov);

			upinit = bns::Vec3F::Normalize(upinit);

			scene.Cameras.push_back(bns::Camera(eyeinit, center, upinit, fov, camera_width, camera_height));
		}
		else if (line.rfind("vertex", 0) == 0)
		{
			const char* ptr_to_vertex = line.c_str();

			bns::I32 end_index = 0;
			bns::Vec3F vert;
			bns::ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.X);
			bns::ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.Y);
			bns::ReadF32FromString(ptr_to_vertex, end_index, &end_index, &vert.Z);
			vertices.push_back(vert);

		}
		else if (line.rfind("pushTransform", 0) == 0)
		{
			if (matrix_stack.empty())
			{
				matrix_stack.push(bns::Mat4x4F::Identity());
			}
			else
			{
				bns::Mat4x4F m = matrix_stack.top();
				m = bns::Mat4x4F::Identity() * m;
				matrix_stack.push(m);
			}
		}
		else if (line.rfind("popTransform", 0) == 0)
		{
			matrix_stack.pop();
		}
		else if (line.rfind("translate", 0) == 0)
		{
			const char* ptr = line.c_str();

			bns::I32 end_index = 0;
			bns::Vec3F v;
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.X);
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.Y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.Z);

			bns::Mat4x4F m = matrix_stack.top();
			matrix_stack.pop();
			m = m * bns::Mat4x4F::Translate(v);
			matrix_stack.push(m);
		}
		else if (line.rfind("scale", 0) == 0)
		{
			const char* ptr = line.c_str();

			bns::I32 end_index = 0;
			bns::Vec3F v;
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.X);
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.Y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.Z);

			bns::Mat4x4F m = matrix_stack.top();
			matrix_stack.pop();
			m = m * bns::Mat4x4F::Scale(v);
			matrix_stack.push(m);
		}
		else if (line.rfind("rotate", 0) == 0)
		{
			const char* ptr = line.c_str();

			bns::I32 end_index = 0;
			bns::Vec3F v;
			bns::F32 theta;
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.X);
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.Y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &v.Z);
			bns::ReadF32FromString(ptr, end_index, &end_index, &theta);

			bns::Mat4x4F m = matrix_stack.top();
			matrix_stack.pop();
			m = m * bns::Mat4x4F::RotationMatrix(bns::Radians(theta), v);
			matrix_stack.push(m);
		}
		else if (line.rfind("ambient", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			bns::ReadF32FromString(ptr, end_index, &end_index, &ambient_color.R);
			bns::ReadF32FromString(ptr, end_index, &end_index, &ambient_color.G);
			bns::ReadF32FromString(ptr, end_index, &end_index, &ambient_color.B);

		}
		else if (line.rfind("emission", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			bns::ReadF32FromString(ptr, end_index, &end_index, &emission_color.R);
			bns::ReadF32FromString(ptr, end_index, &end_index, &emission_color.G);
			bns::ReadF32FromString(ptr, end_index, &end_index, &emission_color.B);
		}
		else if (line.rfind("diffuse", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			bns::ReadF32FromString(ptr, end_index, &end_index, &diffuse_color.R);
			bns::ReadF32FromString(ptr, end_index, &end_index, &diffuse_color.G);
			bns::ReadF32FromString(ptr, end_index, &end_index, &diffuse_color.B);
		}
		else if (line.rfind("specular", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			bns::ReadF32FromString(ptr, end_index, &end_index, &specular_color.R);
			bns::ReadF32FromString(ptr, end_index, &end_index, &specular_color.G);
			bns::ReadF32FromString(ptr, end_index, &end_index, &specular_color.B);
		}
		else if (line.rfind("shininess", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			bns::ReadF32FromString(ptr, end_index, &end_index, &shininess);
		}
		else if (line.rfind("attenuation", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			bns::ReadF32FromString(ptr, end_index, &end_index, &attenuation_constant);
			bns::ReadF32FromString(ptr, end_index, &end_index, &attenuation_linear);
			bns::ReadF32FromString(ptr, end_index, &end_index, &attenuation_quadratic);


		}
		else if (line.rfind("tri", 0) == 0)
		{
			const char* ptr = line.c_str();

			bns::I32 end_index = 0;
			bns::I32 i_1 = 0;
			bns::I32 i_2 = 0;
			bns::I32 i_3 = 0;
			bns::Vec3F vert;
			bns::ReadI32FromString(ptr, end_index, &end_index, &i_1);
			bns::ReadI32FromString(ptr, end_index, &end_index, &i_2);
			bns::ReadI32FromString(ptr, end_index, &end_index, &i_3);

			bns::Vec3F A = vertices[i_1];
			bns::Vec3F B = vertices[i_2];
			bns::Vec3F C = vertices[i_3];

			scene.Shapes.emplace_back(new bns::TriangleShape());
			bns::TriangleShape* ptr_triangle = static_cast<bns::TriangleShape*>(scene.Shapes.back());
			ptr_triangle->Triangle = { A, B, C };
			ptr_triangle->SetTransform(matrix_stack.top());
			ptr_triangle->Material.Ambient = ambient_color;
			ptr_triangle->Material.Emission = emission_color;
			ptr_triangle->Material.Diffuse = diffuse_color;
			ptr_triangle->Material.Specular = specular_color;
			ptr_triangle->Material.Shininess = shininess;
		}
		else if (line.rfind("sphere", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;
			bns::F32 x = 0;
			bns::F32 y = 0;
			bns::F32 z = 0;
			bns::F32 r = 0; // radius
			bns::Vec3F vert;
			bns::ReadF32FromString(ptr, end_index, &end_index, &x);
			bns::ReadF32FromString(ptr, end_index, &end_index, &y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &z);
			bns::ReadF32FromString(ptr, end_index, &end_index, &r);

			scene.Shapes.emplace_back(new bns::SphereShape());
			bns::SphereShape* ptr_sphere = static_cast<bns::SphereShape*>(scene.Shapes.back());
			ptr_sphere->Sphere.Position = bns::Vec3F(x, y, z);
			ptr_sphere->Sphere.Radius = r;
			ptr_sphere->SetTransform(matrix_stack.top());
			ptr_sphere->Material.Ambient = ambient_color;
			ptr_sphere->Material.Diffuse = diffuse_color;
			ptr_sphere->Material.Emission = emission_color;
			ptr_sphere->Material.Specular = specular_color;
			ptr_sphere->Material.Shininess = shininess;
			ptr_sphere->Material.RefractiveIndex = 1.5f;
		}
		else if (line.rfind("cube", 0) == 0)
		{
			scene.Shapes.emplace_back(new bns::CubeShape());
			bns::CubeShape* ptr_cube = static_cast<bns::CubeShape*>(scene.Shapes.back());
			ptr_cube->SetTransform(matrix_stack.top());
			ptr_cube->Material.Ambient = ambient_color;
			ptr_cube->Material.Diffuse = diffuse_color;
			ptr_cube->Material.Emission = emission_color;
			ptr_cube->Material.Specular = specular_color;
			ptr_cube->Material.Shininess = shininess;
		}
		else if (line.rfind("plane", 0) == 0)
		{
			scene.Shapes.emplace_back(new bns::PlaneShape());
			bns::PlaneShape* ptr_plane = static_cast<bns::PlaneShape*>(scene.Shapes.back());
			ptr_plane->SetTransform(matrix_stack.top());
			ptr_plane->Material.Ambient = ambient_color;
			ptr_plane->Material.Diffuse = diffuse_color;
			ptr_plane->Material.Emission = emission_color;
			ptr_plane->Material.Specular = specular_color;
			ptr_plane->Material.Shininess = shininess;
		}
		else if (line.rfind("cylinder", 0) == 0)
		{
			scene.Shapes.emplace_back(new bns::CylinderShape());
			bns::CylinderShape* ptr_plane = static_cast<bns::CylinderShape*>(scene.Shapes.back());
			ptr_plane->SetTransform(matrix_stack.top());
			ptr_plane->Material.Ambient = ambient_color;
			ptr_plane->Material.Diffuse = diffuse_color;
			ptr_plane->Material.Emission = emission_color;
			ptr_plane->Material.Specular = specular_color;
			ptr_plane->Material.Shininess = shininess;
		}
		else if (line.rfind("point", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			scene.Lights.emplace_back(new bns::PointLight());
			bns::PointLight& light = *(static_cast<bns::PointLight*>(scene.Lights.back()));
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Position.X);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Position.Y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Position.Z);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Color.R);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Color.G);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Color.B);

			light.Attenuation.Constant = attenuation_constant;
			light.Attenuation.Linear = attenuation_linear;
			light.Attenuation.Quadratic = attenuation_quadratic;
		}
		else if (line.rfind("directional", 0) == 0)
		{
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;

			scene.Lights.emplace_back(new bns::DirectionalLight());
			bns::DirectionalLight& light = *(static_cast<bns::DirectionalLight*>(scene.Lights.back()));
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Direction.X);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Direction.Y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Direction.Z);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Color.R);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Color.G);
			bns::ReadF32FromString(ptr, end_index, &end_index, &light.Color.B);
		}
	}


	bns::FreeFileContents(file);
}