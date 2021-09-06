#include "bones_loaders.hpp"
#include <vector>
#include <sstream>
#include <stack>
#include "bones_math.hpp"
#include "bones_files.hpp"
#include "bones_char_string.hpp"
#include "bones_scene_common.hpp"
#include <iostream>

/// <summary>
/// Helper method which simply checks if shape should be added to group or scene.
/// Also sets the shape transform matrix from stack.
/// </summary>
void MovePtrToGroupOrScene(bns::RaytracerScene& scene, std::stack<bns::Mat4x4F>& matrix_stack, bns::GroupShape* group_shape, bns::BaseShape* ptr_shape)
{
	if (group_shape != nullptr)
	{
		group_shape->AddChild(ptr_shape);
	}
	else
	{
		scene.Shapes.push_back(ptr_shape);
	}
	ptr_shape->SetTransform(matrix_stack.top());
}

void bns::LoadRaytracerSceneFromBnsFileFormat(const char* filename, bns::RaytracerScene& scene)
{
	std::vector<bns::Vec3F> vertices;

	bns::I32 camera_width = 0;
	bns::I32 camera_height = 0;

	bns::FileContents* file = bns::ReadAndCloseFile(filename);

	std::string line;
	std::stringstream file_stream(file->Contents);

	std::stack<bns::Mat4x4F> matrix_stack;
	matrix_stack.push(bns::Mat4x4F::Identity());

	bns::GroupShape* group_shape = nullptr;

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
		else if (line.rfind("size", 0) == 0)
		{
			bns::I32 index = 0;

			bns::ReadI32FromString(line.c_str(), index, &index, &camera_width);
			bns::ReadI32FromString(line.c_str(), index, &index, &camera_height);
		}
		else if (line.rfind("camera", 0) == 0)
		{
			bns::I32 index = 0;
			const char* ptr_to_camera = line.c_str();

			bns::Vec3F eyeinit;
			bns::ReadF32FromString(ptr_to_camera, index, &index, &eyeinit.X);
			bns::ReadF32FromString(ptr_to_camera, index, &index, &eyeinit.Y);
			bns::ReadF32FromString(ptr_to_camera, index, &index, &eyeinit.Z);

			bns::Vec3F center;
			bns::ReadF32FromString(ptr_to_camera, index, &index, &center.X);
			bns::ReadF32FromString(ptr_to_camera, index, &index, &center.Y);
			bns::ReadF32FromString(ptr_to_camera, index, &index, &center.Z);

			bns::Vec3F upinit;
			bns::ReadF32FromString(ptr_to_camera, index, &index, &upinit.X);
			bns::ReadF32FromString(ptr_to_camera, index, &index, &upinit.Y);
			bns::ReadF32FromString(ptr_to_camera, index, &index, &upinit.Z);

			bns::F32 fov;
			bns::ReadF32FromString(ptr_to_camera, index, &index, &fov);

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

			bns::TriangleShape* ptr_triangle = new bns::TriangleShape(A,B,C);
			MovePtrToGroupOrScene(scene, matrix_stack, group_shape, ptr_triangle);
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


			bns::SphereShape* ptr_sphere = new bns::SphereShape();
			MovePtrToGroupOrScene(scene, matrix_stack, group_shape, ptr_sphere);
			ptr_sphere->Sphere.Position = bns::Vec3F(x, y, z);
			ptr_sphere->Sphere.Radius = r;
			ptr_sphere->Material.Ambient = ambient_color;
			ptr_sphere->Material.Diffuse = diffuse_color;
			ptr_sphere->Material.Emission = emission_color;
			ptr_sphere->Material.Specular = specular_color;
			ptr_sphere->Material.Shininess = shininess;
			ptr_sphere->Material.RefractiveIndex = 1.5f;
		}
		else if (line.rfind("cube", 0) == 0)
		{
			bns::CubeShape* ptr_cube = new bns::CubeShape();
			MovePtrToGroupOrScene(scene, matrix_stack, group_shape, ptr_cube);
			ptr_cube->Material.Ambient = ambient_color;
			ptr_cube->Material.Diffuse = diffuse_color;
			ptr_cube->Material.Emission = emission_color;
			ptr_cube->Material.Specular = specular_color;
			ptr_cube->Material.Shininess = shininess;
		}
		else if (line.rfind("plane", 0) == 0)
		{
			bns::PlaneShape* ptr_plane = new bns::PlaneShape();
			MovePtrToGroupOrScene(scene, matrix_stack, group_shape, ptr_plane);
			ptr_plane->Material.Ambient = ambient_color;
			ptr_plane->Material.Diffuse = diffuse_color;
			ptr_plane->Material.Emission = emission_color;
			ptr_plane->Material.Specular = specular_color;
			ptr_plane->Material.Shininess = shininess;
		}
		else if (line.rfind("cylinder", 0) == 0)
		{
			bns::CylinderShape* ptr_cylinder = new bns::CylinderShape();
			MovePtrToGroupOrScene(scene, matrix_stack, group_shape, ptr_cylinder);
			ptr_cylinder->Material.Ambient = ambient_color;
			ptr_cylinder->Material.Diffuse = diffuse_color;
			ptr_cylinder->Material.Emission = emission_color;
			ptr_cylinder->Material.Specular = specular_color;
			ptr_cylinder->Material.Shininess = shininess;
		}
		else if (line.rfind("pushGroup") == 0)
		{
			bns::GroupShape* new_group = new bns::GroupShape();

			// if group shape already exists, current group is child therefore it goes to group shape.
			if (group_shape != nullptr)
			{
				new_group->Parent = group_shape;
				group_shape->AddChild(new_group);
			}
			// if there is not group shape, group goes directly in shapes.
			else
			{
				scene.Shapes.push_back(new_group);
				group_shape = new_group;
			}
			new_group->SetTransform(matrix_stack.top());
		}
		else if (line.rfind("popGroup") == 0)
		{
			// if there is a group shape
			if (group_shape != nullptr)
			{
				// if there is parent, assign it to be new group_shape
				if (group_shape->Parent != nullptr)
				{
					group_shape = group_shape->Parent;
				}
				// otherwise there is no more grouping, assign pointer to null.
				else
				{
					group_shape = nullptr;
				}
			}
			else
			{
				std::wcout <<
					"Calling 'popGroup' but there is no group to be poped! 'pushGroup' must be called before 'popGroup' if there is any group to be added."
					<< std::endl;
			}
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

void bns::LoadSceneFromObjFileFormat(const char* filename, RaytracerScene& scene)
{
	bns::Camera camera(bns::Vec3F(0.f, 1.5f, 0.f), bns::Vec3F(0.f, 1.f, 0.f), bns::Vec3F(0.f, 1.f, 0.f), 65, 200, 200);
	scene.Cameras.push_back(camera);

	std::vector<bns::Vec3F> vertices;
	std::vector < bns::I32> indices;

	bns::FileContents* file = bns::ReadAndCloseFile(filename);

	std::string line;
	std::stringstream file_stream(file->Contents);

	// TODO: later it should be mesh shape  probably? 
	bns::GroupShape* group_shape = new GroupShape();

	// TODO: default lights 
	bns::PointLight *light = new bns::PointLight();
	light->Color = bns::ColorF(1, 1, 1);
	light->Attenuation.Constant = 1;
	scene.Lights.push_back(light);


	bns::ColorF ambient_color = bns::ColorF(.2f, .2f, .2f);
	//bns::ColorF emission_color;
	//bns::ColorF diffuse_color;
	//bns::ColorF specular_color;
	//bns::F32 shininess = 0.0f;
	//bns::F32 attenuation_constant = 1.0f;
	//bns::F32 attenuation_linear = 0.0f;
	//bns::F32 attenuation_quadratic = 0.0f;

	while (std::getline(file_stream, line))
	{
		if (line.rfind("vt", 0) == 0)
		{
			// TODO: tex coords
		}
		else if (line.rfind("vn", 0) == 0) 
		{
			// TODO: normals
		}
		else if (line.rfind("v", 0) == 0)
		{
			// Vertices
			bns::Vec3F vertex;
			const char* ptr = line.c_str();
			bns::I32 end_index = 0;
			bns::ReadF32FromString(ptr, end_index, &end_index, &vertex.X);
			bns::ReadF32FromString(ptr, end_index, &end_index, &vertex.Y);
			bns::ReadF32FromString(ptr, end_index, &end_index, &vertex.Z);
			vertices.push_back(vertex);
		}
		else if (line.rfind("f", 0) == 0)
		{
			// Faces
			bns::I32 end_index = 0;

			// f v1[/vt1][/vn1] v2[/vt2][/vn2] v3[/vt3][/vn3]
			bns::I32 index = 0;
			bns::I32 i;
			while (bns::ReadI32FromString(line.c_str(), end_index, &end_index, &i))
			{
				if (index == 0)
				{
					indices.push_back(i);
				}

				index++;
				if (index > 1)
				{
					index = 0;
				}
			}

		}
	}

	for (size_t i = 0; i < indices.size(); i += 3)
	{
		bns::Vec3F a = vertices[indices[i]-1];
		bns::Vec3F b = vertices[indices[i + 1]-1];
		bns::Vec3F c = vertices[indices[i + 2]-1];

		bns::TriangleShape*ptr = new bns::TriangleShape(a, b, c);
		//group_shape->AddChild();
		ptr->Material.Specular = bns::ColorF(0.3f, 0.3f, 0.3f);
		ptr->Material.Diffuse = bns::ColorF(0.9, 0.9, 0.9);
		ptr->Material.Ambient = ambient_color;
		ptr->SetTransform(bns::Mat4x4F::Translate(0.0f, 0.0f, -10.0f) * bns::Mat4x4F::Scale(0.1f, 0.1f, 0.1f));
		scene.Shapes.push_back(ptr);
	}

	scene.Shapes.push_back(group_shape);
}
