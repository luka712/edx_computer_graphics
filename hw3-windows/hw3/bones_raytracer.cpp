#include "bones_raytracer.hpp"
#include <iostream>
#include <thread>
#include <vector>
#include <list>

bns::Camera::Camera(bns::Vec3F look_from, bns::Vec3F look_at, bns::Vec3F up, F32 fov, U32 screen_width, U32 screen_height)
	: LookFrom(look_from), LookAt(look_at), Up(up), FOV(fov), ScreenWidth(screen_width), ScreenHeight(screen_height)
{
}

F32 bns::Camera::AspectRatio() const
{
	F32 result = static_cast<F32>(ScreenWidth) / static_cast<F32>(ScreenHeight);
	return result;
}

F32 bns::Camera::FOVInRadians() const
{
	F32 result = bns::Radians(FOV);
	return result;
}

bns::RayF bns::RayThroughPixel(Camera cam, I32 pixel_x, I32 pixel_y)
{
	// Steps
	// 1. Create a coordinate frame for the camera
	// 2. Define a rotation matrix
	// 3. Apply appropriate translation for camera ( eye ) location

	//      a          b x w
	// w = ---    u = -------       v = w x u
	//    ||a||     || b x w ||

	// a = eye - center
	bns::Vec3F eye = cam.LookFrom;
	bns::Vec3F center = cam.LookAt;
	bns::Vec3F Up = cam.Up;

	bns::Vec3F a = eye - center;
	bns::Vec3F w = bns::Vec3F::Normalize(a);

	bns::Vec3F b = bns::Vec3F::Normalize(Up);

	bns::Vec3F b_cross_w = bns::Vec3F::Cross(b, w);
	bns::Vec3F b_cross_w_unit = bns::Vec3F::Normalize(b_cross_w);

	bns::Vec3F u = b_cross_w_unit;

	bns::Vec3F v = bns::Vec3F::Cross(w, u);

	F32 half_w = cam.ScreenWidth / 2.0f;
	F32 half_h = cam.ScreenHeight / 2.0f;

	F32 _i = pixel_x + 0.5f;
	F32 _j = pixel_y + 0.5f;

	F32 fov = cam.FOVInRadians();

	// https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	F32 half_view = bns::Tan(fov / 2.0f);

	F32 aspect_ratio = cam.AspectRatio();

	// TODO: fix for aspect ration
	// in this assignemnt x is wider, therefore fix alpha, otherwise beta is to be fixed.

	F32 alpha = half_view * aspect_ratio * ((_i - half_w) / half_w);
	F32 beta = half_view * ((half_h - _j) / half_h);

	bns::Vec3F dir = u * alpha + v * beta - w;

	dir.Normalize();

	// TODO: vec3 to point3f
	bns::RayF ray(bns::Point3F(eye.X, eye.Y, eye.Z), dir);
	return ray;
}

bns::Intersections bns::GetIntersections(const bns::RayF& ray, bns::BaseShape** Shapes, U32 Count)
{
	bns::Intersections result;

	F32 min_dist = MAX_F32;
	const BaseShape* ptr_hit_shape = nullptr;
	for (U32 i = 0; i < Count; i++)
	{
		const BaseShape* shape = Shapes[i];

		bns::RayF transformed_ray = ray * shape->GetInverseTransform();

		// Note that intersections are passed and filled in intersection distance
		shape->IntersectRayAndObjects(transformed_ray, result);
	}

	return result;
}

bns::ColorF** bns::AllocateColors(const bns::Camera& camera)
{
	bns::ColorF** pixels = (bns::ColorF**)malloc(camera.ScreenHeight * sizeof(bns::ColorF*));
	for (U32 i = 0; i < camera.ScreenHeight; i++)
	{
		pixels[i] = (bns::ColorF*)malloc(camera.ScreenWidth * sizeof(bns::ColorF));
	}
	return pixels;
}

bns::ColorF bns::ColorAt(const bns::RayF& ray,
	bns::BaseShape** shapes, U32 count_of_shapes,
	bns::BaseLight** lights, U32 count_of_lights,
	I32 remaining)
{
	if (remaining == 0)
	{
		return bns::ColorF::Black();
	}

	// 2. Intersect that ray with all objects 
	bns::Intersections intersections = GetIntersections(ray, shapes, count_of_shapes);

	// 3. Get intersections that's hit by ray 
	bns::Intersection* intersection = intersections.Hit();

	if (intersection != nullptr)
	{
		bns::Computations comp = bns::PrepareComputations(*intersection, ray, intersections);

		// B) is done here. Shadow rays are fired in lighting model as well as Illumination model
		bns::ColorF surface_color = BlinnPhongReflectionModel(comp, shapes, count_of_shapes, lights, count_of_lights);

		// C) is done here. Secondary rays are fired off to trace reflected rays.
		bns::ColorF reflected_color = ReflectedColor(comp, shapes, count_of_shapes, lights, count_of_lights, remaining);

		// C) is done here. Secondary rays are fired off to trace refracted rays.
		bns::ColorF refracted_color = RefractedColor(comp, shapes, count_of_shapes, lights, count_of_lights, remaining);

		return surface_color + reflected_color + refracted_color;
	}

	return bns::ColorF::Black();
}

// Local function
void RayTraceFromToHorizontalLine(U32 start_from_line, U32 end_on_line, const bns::Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::BaseLight** lights, U32 count_of_lights, bns::ColorF** colors_to_fill)
{
	for (U32 j = start_from_line; j < end_on_line; j++)
	{
		F32 percent = static_cast<F32>(j - start_from_line) / static_cast<F32>(end_on_line - start_from_line) * 100.0f;

		std::cout << static_cast<I32>(percent) << "/" << 100 << std::endl;
		for (size_t i = 0; i < cam.ScreenWidth; i++)
		{
			//std::cout << i << "/" << cam.ScreenWidth << std::endl;

			// recursive ray tracing 
			// A: Trace primary eye ray, find intersection
			// B: Trace secondary shadow rays to all lights. Color = Visible ? Illumination Model : 0
			// C: trace reflected ray: Color += reflectivity * Color of reflected ray


			// 1. The primary ray. 
			bns::RayF primary_ray = bns::RayThroughPixel(cam, i, j);

			bns::ColorF color = ColorAt(primary_ray, shapes, count_of_shapes, lights, count_of_lights, 5);
			color.A = 1.0f;

			// B and C step inside color at.
			colors_to_fill[j][i] = color;
		}
	}
}

void bns::ThreadedRayTrace(const Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::BaseLight** lights, U32 count_of_lights, bns::ColorF** colors_to_fill)
{
	std::thread threads[10];

	U32 increment = cam.ScreenHeight / 10;
	U32 thread_index = 0;
	for (U32 i = 0; i < cam.ScreenHeight; i += increment)
	{
		threads[thread_index++] = std::thread(RayTraceFromToHorizontalLine, i, i + increment, cam, shapes, count_of_shapes, lights, count_of_lights, colors_to_fill);
	}

	for (size_t i = 0; i < 10; i++)
	{
		threads[i].join();
	}
}

void bns::RayTrace(const Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::BaseLight** lights, U32 count_of_lights, bns::ColorF** colors_to_fill)
{
	for (size_t j = 0; j < cam.ScreenHeight; j++)
	{
		std::cout << j << "/" << cam.ScreenHeight << std::endl;
		for (size_t i = 0; i < cam.ScreenWidth; i++)
		{
			//std::cout << i << "/" << cam.ScreenWidth << std::endl;

			// recursive ray tracing 
			// A: Trace primary eye ray, find intersection
			// B: Trace secondary shadow rays to all lights. Color = Visible ? Illumination Model : 0
			// C: trace reflected ray: Color += reflectivity * Color of reflected ray


			// 1. The primary ray. 
			bns::RayF primary_ray = bns::RayThroughPixel(cam, i, j);

			bns::ColorF color = ColorAt(primary_ray, shapes, count_of_shapes, lights, count_of_lights, 5);
			color.A = 1.0f;

			// B and C step inside color at.
			colors_to_fill[j][i] = color;
		}
	}
}

bns::Computations bns::PrepareComputations(const Intersection& intersection, const RayF& ray, const Intersections& other)
{
	bns::Computations result;

	result.T = intersection.TMinDist;
	result.Shape = intersection.HitShape;
	result.RayOrigin = ray.Origin;
	result.RayDirection = ray.Direction.ToVec3F();
	result.WorldPoint = ray.Origin + ray.Direction * intersection.TMinDist;
	result.LocalPoint = intersection.HitShape->GetInverseTransform() * result.WorldPoint;
	result.Eye = (ray.Direction * -1.0f).ToVec3F();
	result.LocalNormal = result.Shape->GetNormalAt(result.LocalPoint.ToPoint3F());
	result.LocalNormal.Normalize();
	result.WorldNormal = (bns::Mat4x4F::Transpose(result.Shape->GetInverseTransform()) * result.LocalNormal.ToVec4F(1.0f)).ToVec3F();
	result.WorldNormal.Normalize();

	if (result.WorldNormal.Dot(result.Eye) < 0)
	{
		result.WorldNormal = result.WorldNormal * -1.0f;
	}

	result.ReflectedVector = bns::Reflect(ray.Direction.ToVec3F(), result.WorldNormal);

	U32 containers_index = 0;
	// Replace with custom array.
	std::list<const BaseShape*> container;

	for (U32 i = 0; i < other.CurrentIntersectionCount; i++)
	{
		const Intersection& _i = other.IntersectionsArray[i];

		// 1. If the intersection it hit, set n1 to refractive index of last object in the containers list.
		// If that list is emptry, then tere is no containing object, and n1 should be set to 1.

		// compare the pointer address
		if (&_i == &intersection)
		{
			// empty, just push index of 1.0
			if (container.empty())
			{
				result.N1 = 1.0;
			}
			else
			{
				// else use refractive index of material.
				result.N1 = container.back()->Material.RefractiveIndex;
			}
		}

		// 2. If the interesections's object is already in the container list, then this intersection must be exiting the object.
		// Remove the objects from container in this case. Otherwise, intersection is entering the object, and the object should be added to the end of list.

		// if container includes intersection object, remove it.
		bool is_included = false;
		for (auto it = container.begin(); it != container.end(); ++it)
		{
			// found object by address
			if (*it == _i.HitShape)
			{
				container.remove(*it);
				is_included = true;
				break;
			}
		}
		if (!is_included)
		{
			container.emplace_back(_i.HitShape);
		}

		// 3. If the intersection is the hit ( _i is intersection ), set n2 to the refractive index of the last object in the containers list.
		// If the list is empty, then again, there is no containing object and n2 should be set to 1.

		// compare the pointer address
		if (&_i == &intersection)
		{
			if (container.empty())
			{
				result.N2 = 1.0f;
			}
			else
			{
				result.N2 = container.back()->Material.RefractiveIndex;
			}

			// 4. If the intersection is hit, terminate the loop here
			break;
		}
	}

	return result;
}

bool bns::IsShadowed(const Computations& comp, const BaseLight& light_param, bns::BaseShape** shapes, U32 count_of_shapes)
{
	if (light_param.Type == bns::LightType::Point)
	{
		const bns::PointLight& light = static_cast<const bns::PointLight&>(light_param);

		bns::Vec3F point_to_light = light.Position - comp.WorldPoint.ToVec3F();
		point_to_light.Normalize();

		// Move origin slightly towards the direction, in order to not intersect ray with self.
		bns::Point4F origin = comp.WorldPoint + (point_to_light * EPSILON).ToPoint4F(0.0f);

		bns::RayF ray(origin, point_to_light);

		// Get all the intersections of ray with shapes.
		bns::Intersections intersections = GetIntersections(ray, shapes, count_of_shapes);

	
		// if there is a hit, object cannot receive light, as there is another object between it and sun
		// also determine if hit objects is behind light
		Intersection* hit = intersections.Hit();
		if (hit)
		{
			// Move origin slightly towards the direction, in order to not intersect ray with self.
			bns::Point4F trace_light_to_hit_origin = light.Position.ToPoint4F(1.0f) + (point_to_light * EPSILON).ToPoint4F(0.0f);

			bns::RayF ray_from_light(trace_light_to_hit_origin, point_to_light);

			// Get all the intersections of ray with shapes. 
			bns::Intersections intersections_ray_from_light = GetIntersections(ray_from_light, shapes, count_of_shapes);

			// Here is important to determine if the hit object is directly behind the light. If there is no hit from light in direction of normal towards light
			// it is fine, it's also fine if some other object is hit as that means that previously hit object is betweel intersection point and light.
		
			Intersection* hit_from_light = intersections_ray_from_light.Hit();
			if (!hit_from_light || hit_from_light != hit)
			{
				return true;
			}
		}		
	}
	else if (light_param.Type == bns::LightType::Directional)
	{
		const bns::DirectionalLight& light = static_cast<const bns::DirectionalLight&>(light_param);

		bns::Vec3F point_to_light = light.Direction;
		point_to_light.Normalize();

		// Move origin slightly towards the direction, in order to not intersect ray with self.
		bns::Point4F origin = comp.WorldPoint + (point_to_light * EPSILON).ToPoint4F(0.0f);

		bns::RayF ray(origin, point_to_light);

		// Get all the intersections of ray with shapes.
		bns::Intersections intersections = GetIntersections(ray, shapes, count_of_shapes);

		// if there is a hit, object cannot receive light, as there is another object between it and sun
		Intersection* hit = intersections.Hit();
		if (hit)
		{
			return true;
		}
	}

	return false;
}

bns::ColorF bns::BlinnPhongReflectionModel(const Computations& computation, bns::BaseShape** shapes, U32 count_of_shapes, BaseLight** lights, U32 count_of_lights)
{
	// One can assume that there is hit here, as it is only called by code if there is a hit.

	const bns::BaseShape* shape = computation.Shape;

	bns::ColorF ambient = shape->Material.Ambient;
	bns::ColorF emissive = shape->Material.Emission;
	bns::ColorF diffuse;
	bns::ColorF specular;

	for (size_t i = 0; i < count_of_lights; i++)
	{
		bns::Vec3F light_dir;
		bns::Vec3F world_normal;

		if (lights[i]->Type == bns::LightType::Point)
		{
			bns::PointLight light = *static_cast<bns::PointLight*>(lights[i]);

			bool is_shadowed = IsShadowed(computation, light, shapes, count_of_shapes);

			if (is_shadowed)
			{
				continue;
			}

			world_normal = computation.WorldNormal;

			light_dir = light.Position - computation.WorldPoint.ToVec3F();
			light_dir.Normalize();

			F32 n_dot_l = world_normal.Dot(light_dir);

			F32 distance = Distance(computation.WorldPoint.ToVec3F(), light.Position);
			F32 attenuation = 1.0f / (light.Attenuation.Constant + light.Attenuation.Linear * distance + light.Attenuation.Quadratic * Pow(distance));

			if (n_dot_l > 0.0f)
			{
				diffuse += shape->Material.Diffuse * light.Color * attenuation * n_dot_l;
			}
		}
		else if (lights[i]->Type == bns::LightType::Directional)
		{
			bns::DirectionalLight light = *static_cast<bns::DirectionalLight*>(lights[i]);

			bool is_shadowed = IsShadowed(computation, light, shapes, count_of_shapes);

			if (is_shadowed)
			{
				continue;
			}

			bns::Vec3F world_normal = computation.WorldNormal;

			light_dir = light.Direction;
			light_dir.Normalize();

			F32 n_dot_l = world_normal.Dot(light_dir);

			if (n_dot_l > 0.0f)
			{
				diffuse += shape->Material.Diffuse * light.Color * n_dot_l;
			}
		}

		// Specular blinn-phong
		bns::Vec3F half_v = bns::Normalize(computation.Eye + light_dir);
		F32 n_dot_h = bns::Dot(world_normal, half_v);

		if (n_dot_h > 0.0f)
		{
			specular += shape->Material.Specular * lights[i]->Color * Pow(n_dot_h, shape->Material.Shininess);
		}

	}

	bns::ColorF result = ambient + emissive + diffuse + specular;
	result.A = 1.0f;
	return result;
}

bns::ColorF bns::ReflectedColor(const Computations& comp,
	bns::BaseShape** shapes, U32 count_of_shapes,
	bns::BaseLight** lights, U32 count_of_lights,
	I32 remaining)
{
	if (remaining <= 0)
	{
		return bns::ColorF::Black();
	}

	// C) shooting of reflected rays.

	// Move by epsilon a bit in direction of normal
	bns::Point4F origin = comp.WorldPoint + comp.ReflectedVector.ToVec4F() * EPSILON;

	// Shoot reflection from a point in direction of normal.
	bns::RayF reflected_ray(origin, comp.ReflectedVector);
	bns::ColorF color = bns::ColorAt(reflected_ray, shapes, count_of_shapes, lights, count_of_lights, --remaining);


	return color * comp.Shape->Material.Specular;
}

bns::ColorF bns::RefractedColor(const Computations& comp,
	bns::BaseShape** shapes, U32 count_of_shapes,
	bns::BaseLight** lights, U32 count_of_lights,
	I32 remaining)
{
	if (remaining < 0.0f)
	{
		return bns::ColorF::Black();
	}

	// Find the ratio of first index of refraction to the second.
	// Inverted from the definition of Snell's Law
	F32 n_ratio = comp.N1 / comp.N2;

	// cos(theta_i) is the same as dot product of two vectors
	F32 cos_i = comp.Eye.Dot(comp.LocalNormal);

	// Find sin(theta_t)^2 via trigonometric identity
	F32 sin2_t = (n_ratio * n_ratio) * (1 - cos_i * cos_i);

	// Find cos(theta_t) via trigonometric identity.
	F32 cos_t = bns::Sqrt(1.0 - sin2_t);

	bns::Vec3F direction_of_refracted_ray = (comp.WorldNormal * (n_ratio * cos_i - cos_t)) - (comp.Eye * (n_ratio));

	bns::RayF refracted_ray = bns::RayF(comp.WorldPoint + direction_of_refracted_ray.ToVec4F() * EPSILON, direction_of_refracted_ray);

	// Multiply by transparency to account for opacity.
	bns::ColorF color = ColorAt(refracted_ray, shapes, count_of_shapes, lights, count_of_lights, --remaining);

	//return this.colorAt(refracted_ray, remaining - 1).multiplyByScalar(c.object.material.transparency);
	return color * comp.Shape->Material.Diffuse;
}

void* bns::ColorsToABGR8888Pixels(const Camera& camera, bns::ColorF** colors)
{
	I32* result = (I32*)malloc(camera.ScreenWidth * camera.ScreenHeight * sizeof(I32));

	U32 pixel_index = 0;
	for (U32 y = 0; y < camera.ScreenHeight; y++)
	{
		for (U32 x = 0; x < camera.ScreenWidth; x++)
		{
			bns::ColorF color = colors[y][x];
			result[pixel_index] = color.ToABGR8888();
			pixel_index++;
		}
	}

	return result;
}

void* bns::ColorsToARGB8888Pixels(const Camera& camera, bns::ColorF** colors)
{
	I32* result = (I32*)malloc(camera.ScreenWidth * camera.ScreenHeight * sizeof(I32));

	U32 pixel_index = 0;
	for (U32 y = 0; y < camera.ScreenHeight; y++)
	{
		for (U32 x = 0; x < camera.ScreenWidth; x++)
		{
			bns::ColorF color = colors[y][x];
			result[pixel_index] = color.ToARGB8888();
			pixel_index++;
		}
	}

	return result;
}



void bns::FreeColors(const bns::Camera& camera, bns::ColorF** colors)
{
	for (U32 i = 0; i < camera.ScreenHeight; i++)
	{
		free(colors[i]);
	}
	free(colors);
}

void bns::FreePixels(void* pixels)
{
	free(pixels);
}


bns::BaseShape::BaseShape()
{
	Material = {};
	Transform = bns::Mat4x4F::Identity();
}

const bns::Mat4x4F& bns::BaseShape::GetTransform() const
{
	return Transform;
}

const bns::Mat4x4F& bns::BaseShape::GetInverseTransform() const
{
	return this->InverseTransform;
}

void bns::BaseShape::SetTransform(bns::Mat4x4F m)
{
	Transform = m;
	InverseTransform = bns::Mat4x4F::Inverse(m);
}


bns::TriangleShape::TriangleShape()
{
	Triangle = bns::TriangleF();
}

void bns::TriangleShape::IntersectRayAndObjects(const bns::RayF& ray, Intersections& fill_intersections) const
{
	F32 result = ray.IntersectionDistanceWithTriangle(this->Triangle);

	if (result != MAX_F32)
	{
		Intersection i;
		i.TMinDist = result;
		i.HitShape = this;
		i.Type = ShapeType::Triangle;
		fill_intersections.Add(i);
	}
}

bns::Vec3F bns::TriangleShape::GetNormalAt(bns::Point3F point) const
{
	bns::Vec3F result = this->Triangle.GetNormal();
	return result;
}

void bns::SphereShape::IntersectRayAndObjects(const bns::RayF& ray, Intersections& fill_intersections) const
{
	F32 t1, t2;
	ray.IntersectionDistanceWithSphere(this->Sphere, &t1, &t2);

	if (t1 != MAX_F32)
	{
		Intersection i1;
		i1.TMinDist = t1;
		i1.HitShape = this;
		i1.Type = ShapeType::Triangle;
		fill_intersections.Add(i1);
	}

	if (t1 != MAX_F32)
	{
		Intersection i2;
		i2.TMinDist = t2;
		i2.HitShape = this;
		i2.Type = ShapeType::Triangle;
		fill_intersections.Add(i2);
	}
}

bns::Vec3F bns::SphereShape::GetNormalAt(bns::Point3F point) const
{
	return point.ToVec3F();
}

bns::SphereShape::SphereShape()
{
	Sphere = bns::SphereF();
}

bns::Intersections::Intersections()
{
	CurrentIntersectionCount = 0;
}

void bns::Intersections::Add(Intersection intersection)
{
	// TODO: move out of here. this will fail at some point, do increase count of intersections if neccessary.
	IntersectionsArray[CurrentIntersectionCount].HitShape = intersection.HitShape;
	IntersectionsArray[CurrentIntersectionCount].TMinDist = intersection.TMinDist;
	IntersectionsArray[CurrentIntersectionCount].Type = intersection.Type;
	CurrentIntersectionCount++;
}

bns::Intersection* bns::Intersections::Hit()
{
	I32 min_hit_index = -1;
	for (I32 i = 0; i < static_cast<I32>(CurrentIntersectionCount); i++)
	{
		if (this->IntersectionsArray[i].TMinDist > 0.0f)
		{
			if (min_hit_index == -1)
			{
				min_hit_index = i;
			}

			if (this->IntersectionsArray[i].TMinDist < this->IntersectionsArray[min_hit_index].TMinDist)
			{
				min_hit_index = i;
			}
		}
	}

	if (min_hit_index > -1)
	{
		return &this->IntersectionsArray[min_hit_index];
	}

	return nullptr;
}

bns::Intersection bns::Intersection::operator=(const Intersection& other)
{
	Intersection i;
	i.HitShape = other.HitShape;
	i.TMinDist = other.TMinDist;
	i.Type = other.Type;
	return i;
}
