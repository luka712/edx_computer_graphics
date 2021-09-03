#include "bones_raytracer.hpp"
#include <iostream>
#include <thread>
#include <vector>
#include <list>

bns::Camera::Camera(bns::Vec3F look_from, bns::Vec3F look_at, bns::Vec3F up, bns::F32 fov, bns::U32 screen_width, bns::U32 screen_height)
	: LookFrom(look_from), LookAt(look_at), Up(up), FOV(fov), ScreenWidth(screen_width), ScreenHeight(screen_height)
{
}

bns::F32 bns::Camera::AspectRatio() const
{
	bns::F32 result = static_cast<bns::F32>(ScreenWidth) / static_cast<bns::F32>(ScreenHeight);
	return result;
}

bns::F32 bns::Camera::FOVInRadians() const
{
	bns::F32 result = bns::Radians(FOV);
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

	bns::F32 half_w = cam.ScreenWidth / 2.0f;
	bns::F32 half_h = cam.ScreenHeight / 2.0f;

	bns::F32 _i = pixel_x + 0.5f;
	bns::F32 _j = pixel_y + 0.5f;

	bns::F32 fov = cam.FOVInRadians();

	// https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	bns::F32 half_view = bns::Tan(fov / 2.0f);

	bns::F32 aspect_ratio = cam.AspectRatio();

	// TODO: fix for aspect ration
	// in this assignemnt x is wider, therefore fix alpha, otherwise beta is to be fixed.

	bns::F32 alpha = half_view * aspect_ratio * ((_i - half_w) / half_w);
	bns::F32 beta = half_view * ((half_h - _j) / half_h);

	bns::Vec3F dir = u * alpha + v * beta - w;

	dir.Normalize();

	// TODO: vec3 to point3f
	bns::RayF ray(bns::Point3F(eye.X, eye.Y, eye.Z), dir);
	return ray;
}

bns::Intersections bns::GetIntersections(const bns::RayF& ray, bns::BaseShape** Shapes, bns::U32 Count)
{
	bns::Intersections result;

	bns::F32 min_dist = bns::MAX_F32;
	const BaseShape* ptr_hit_shape = nullptr;
	for (bns::U32 i = 0; i < Count; i++)
	{
		const BaseShape* shape = Shapes[i];

		bns::RayF transformed_ray = ray * shape->GetInverseTransform();

		// Note that intersections are passed and filled in intersection distance
		shape->IntersectWithRay(transformed_ray, result);
	}

	return result;
}

bns::ColorF** bns::AllocateColors(const bns::Camera& camera)
{
	bns::ColorF** pixels = (bns::ColorF**)malloc(camera.ScreenHeight * sizeof(bns::ColorF*));
	for (bns::U32 i = 0; i < camera.ScreenHeight; i++)
	{
		pixels[i] = (bns::ColorF*)malloc(camera.ScreenWidth * sizeof(bns::ColorF));
	}
	return pixels;
}

bns::ColorF bns::ColorAt(const bns::RayF& ray,
	bns::BaseShape** shapes, bns::U32 count_of_shapes,
	bns::BaseLight** lights, bns::U32 count_of_lights,
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
void RayTraceFromToHorizontalLine(bns::U32 start_from_line, bns::U32 end_on_line, const bns::Camera& cam, bns::BaseShape** shapes, bns::U32 count_of_shapes, bns::BaseLight** lights, bns::U32 count_of_lights, bns::ColorF** colors_to_fill)
{
	for (bns::U32 j = start_from_line; j < end_on_line; j++)
	{
		bns::F32 percent = static_cast<bns::F32>(j - start_from_line) / static_cast<bns::F32>(end_on_line - start_from_line) * 100.0f;

		std::cout << static_cast<bns::I32>(percent) << "/" << 100 << std::endl;
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

void bns::ThreadedRayTrace(const Camera& cam, bns::BaseShape** shapes, bns::U32 count_of_shapes, bns::BaseLight** lights, bns::U32 count_of_lights, bns::ColorF** colors_to_fill)
{
	std::thread threads[10];

	bns::U32 increment = cam.ScreenHeight / 10;
	bns::U32 thread_index = 0;
	for (bns::U32 i = 0; i < cam.ScreenHeight; i += increment)
	{
		threads[thread_index++] = std::thread(RayTraceFromToHorizontalLine, i, i + increment, cam, shapes, count_of_shapes, lights, count_of_lights, colors_to_fill);
	}

	for (size_t i = 0; i < 10; i++)
	{
		threads[i].join();
	}
}

void bns::RayTrace(const Camera& cam, bns::BaseShape** shapes, bns::U32 count_of_shapes, bns::BaseLight** lights, bns::U32 count_of_lights, bns::ColorF** colors_to_fill)
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
	result.WorldPoint = (ray.Origin + ray.Direction * intersection.TMinDist);
	result.LocalPoint = WorldPointToLocalPoint(*result.Shape, result.WorldPoint);
	result.Eye = (ray.Direction * -1.0f).ToVec3F();
	result.LocalNormal = result.Shape->GetLocalNormalAt(result.LocalPoint.ToPoint3F());
	result.LocalNormal.Normalize();
	result.WorldNormal = GetWorldNormalAt(*result.Shape, result.LocalNormal);
	result.WorldNormal.Normalize();

	if (result.WorldNormal.Dot(result.Eye) < 0)
	{
		result.WorldNormal = result.WorldNormal * -1.0f;
	}

	result.ReflectedVector = bns::Reflect(ray.Direction.ToVec3F(), result.WorldNormal);

	bns::U32 containers_index = 0;
	// Replace with custom array.
	std::list<const BaseShape*> container;

	for (bns::U32 i = 0; i < other.CurrentIntersectionCount; i++)
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

bool bns::IsShadowed(const Computations& comp, const BaseLight& light_param, bns::BaseShape** shapes, bns::U32 count_of_shapes)
{
	if (light_param.Type == bns::LightType::Point)
	{
		const bns::PointLight& light = static_cast<const bns::PointLight&>(light_param);

		bns::Vec3F point_to_light = light.Position - comp.WorldPoint.ToVec3F();
		F32 distance = point_to_light.Length();
		point_to_light.Normalize();

		// Move origin slightly towards the direction, in order to not intersect ray with self.
		bns::Point4F origin = comp.WorldPoint + (point_to_light * EPSILON).ToPoint4F(0.0f);

		bns::RayF ray(origin, point_to_light);

		// Get all the intersections of ray with shapes.
		bns::Intersections intersections = GetIntersections(ray, shapes, count_of_shapes);


		// if there is a hit, object cannot receive light, as there is another object between it and sun
		// also determine if hit objects is behind light
		Intersection* hit = intersections.Hit();
		if (hit && hit->TMinDist < distance)
		{
			return true;


			// Move origin slightly towards the direction, in order to not intersect ray with self.
			bns::Point4F trace_light_to_hit_origin = light.Position.ToPoint4F(1.0f) + (point_to_light * EPSILON).ToPoint4F(0.0f);

			bns::RayF ray_from_light(trace_light_to_hit_origin, point_to_light);

			// Get all the intersections of ray with shapes. 
			bns::Intersections intersections_ray_from_light = GetIntersections(ray_from_light, shapes, count_of_shapes);

			// Here is important to determine if the hit object is directly behind the light. If there is no hit from light in direction of normal towards light
			// it is fine, it's also fine if some other object is hit as that means that previously hit object is betweel intersection point and light.

			Intersection* hit_from_light = intersections_ray_from_light.Hit();
			if (!hit_from_light && hit_from_light != hit)
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

bns::ColorF bns::BlinnPhongReflectionModel(const Computations& computation, bns::BaseShape** shapes, bns::U32 count_of_shapes, BaseLight** lights, bns::U32 count_of_lights)
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

			bns::F32 n_dot_l = world_normal.Dot(light_dir);

			bns::F32 distance = Distance(computation.WorldPoint.ToVec3F(), light.Position);
			bns::F32 attenuation = 1.0f / (light.Attenuation.Constant + light.Attenuation.Linear * distance + light.Attenuation.Quadratic * Pow(distance));

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

			bns::F32 n_dot_l = world_normal.Dot(light_dir);

			if (n_dot_l > 0.0f)
			{
				diffuse += shape->Material.Diffuse * light.Color * n_dot_l;
			}
		}

		// Specular blinn-phong
		bns::Vec3F half_v = bns::Normalize(computation.Eye + light_dir);
		bns::F32 n_dot_h = bns::Dot(world_normal, half_v);

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
	bns::BaseShape** shapes, bns::U32 count_of_shapes,
	bns::BaseLight** lights, bns::U32 count_of_lights,
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
	bns::BaseShape** shapes, bns::U32 count_of_shapes,
	bns::BaseLight** lights, bns::U32 count_of_lights,
	I32 remaining)
{
	if (remaining < 0.0f)
	{
		return bns::ColorF::Black();
	}

	// Find the ratio of first index of refraction to the second.
	// Inverted from the definition of Snell's Law
	bns::F32 n_ratio = comp.N1 / comp.N2;

	// cos(theta_i) is the same as dot product of two vectors
	bns::F32 cos_i = comp.Eye.Dot(comp.LocalNormal);

	// Find sin(theta_t)^2 via trigonometric identity
	bns::F32 sin2_t = (n_ratio * n_ratio) * (1 - cos_i * cos_i);

	// Find cos(theta_t) via trigonometric identity.
	bns::F32 cos_t = bns::Sqrt(1.0 - sin2_t);

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

	bns::U32 pixel_index = 0;
	for (bns::U32 y = 0; y < camera.ScreenHeight; y++)
	{
		for (bns::U32 x = 0; x < camera.ScreenWidth; x++)
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

	bns::U32 pixel_index = 0;
	for (bns::U32 y = 0; y < camera.ScreenHeight; y++)
	{
		for (bns::U32 x = 0; x < camera.ScreenWidth; x++)
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
	for (bns::U32 i = 0; i < camera.ScreenHeight; i++)
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
	Parent = nullptr;
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

void bns::TriangleShape::IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const
{
	bns::F32 result = ray.IntersectionDistanceWithTriangle(this->Triangle);

	if (result != bns::MAX_F32)
	{
		Intersection i;
		i.TMinDist = result;
		i.HitShape = this;
		i.Type = ShapeType::Triangle;
		fill_intersections.Add(i);
	}
}

bns::Vec3F bns::TriangleShape::GetLocalNormalAt(bns::Point3F point) const
{
	bns::Vec3F result = this->Triangle.GetNormal();
	return result;
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

#pragma region SPHERE SHAPE  

void bns::SphereShape::IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const
{
	bns::F32 t1, t2;
	ray.IntersectionDistanceWithSphere(this->Sphere, &t1, &t2);

	if (t1 != bns::INFINITY_F32)
	{
		Intersection i1;
		i1.TMinDist = t1;
		i1.HitShape = this;
		i1.Type = ShapeType::Triangle;
		fill_intersections.Add(i1);
	}

	if (t2 != bns::INFINITY_F32)
	{
		Intersection i2;
		i2.TMinDist = t2;
		i2.HitShape = this;
		i2.Type = ShapeType::Triangle;
		fill_intersections.Add(i2);
	}
}

bns::Vec3F bns::SphereShape::GetLocalNormalAt(bns::Point3F point) const
{
	return point.ToVec3F();
}

bns::SphereShape::SphereShape()
{
	Sphere = bns::SphereF();
}

#pragma endregion

#pragma region PLANE SHAPE

void bns::PlaneShape::IntersectWithRay(const RayF& ray, Intersections& fill_intersections) const
{
	if (Abs(ray.Direction.Y) < EPSILON)
	{
		return;
	}

	F32 t = -ray.Origin.Y / ray.Direction.Y;
	
	Intersection i;
	i.TMinDist = t;
	i.HitShape = this;
	i.Type = ShapeType::Plane;
	fill_intersections.Add(i);
}

bns::Vec3F bns::PlaneShape::GetLocalNormalAt(bns::Point3F point) const
{
	return Vec3F(0.0f, 1.0f, 0.0f);
}

#pragma endregion

#pragma region CUBE SHAPE

bns::CubeShape::CubeShape()
{
}

void bns::CubeShape::CheckAxis(F32 origin, F32 direction, F32* tmin, F32* tmax) const
{
	F32 tmin_numerator = (-1.0f - origin);
	F32 tmax_numerator = (1.0f - origin);


	if (Abs(direction) >= EPSILON)
	{
		*tmin = tmin_numerator / direction;
		*tmax = tmax_numerator / direction;
	}
	else
	{
		*tmin = tmin_numerator * bns::INFINITY_F32;
		*tmax = tmax_numerator * bns::INFINITY_F32;
	}


	if (*tmin > *tmax)
	{
		bns::Swap(tmin, tmax);
	}
}

bns::Vec3F bns::CubeShape::GetLocalNormalAt(bns::Point3F point) const
{
	F32 abs_x = Abs(point.X);
	F32 abs_y = Abs(point.Y);
	F32 abs_z = Abs(point.Z);

	F32 max = Max(abs_x, Max(abs_y, abs_z));

	bns::Vec3F result = bns::Vec3F(0, 0, point.Z);
	if (max == abs_x)
	{
		result = bns::Vec3F(point.X, 0.0f, 0.0f);
	}
	else if (max == abs_y)
	{
		result = bns::Vec3F(0.0f, point.Y, 0.0f);
	}

	return result;
}

void bns::CubeShape::IntersectWithRay(const bns::RayF& ray, Intersections& fill_intersections) const
{
	F32 xtmin, xtmax, ytmin, ytmax, ztmin, ztmax;
	CheckAxis(ray.Origin.X, ray.Direction.X, &xtmin, &xtmax);
	CheckAxis(ray.Origin.Y, ray.Direction.Y, &ytmin, &ytmax);
	CheckAxis(ray.Origin.Z, ray.Direction.Z, &ztmin, &ztmax);

	F32 tmin = bns::Max(xtmin, bns::Max(ytmin, ztmin));
	F32 tmax = bns::Min(xtmax, bns::Min(ytmax, ztmax));

	if (tmax >= tmin)
	{
		if (tmin != bns::INFINITY_F32 && tmin != bns::NEGATIVE_INFINITY_F32)
		{
			bns::Intersection i1;
			i1.HitShape = this;
			i1.TMinDist = tmin;
			i1.Type = ShapeType::Cube;
			fill_intersections.Add(i1);
		}

		if (tmax != bns::INFINITY_F32 && tmax != bns::NEGATIVE_INFINITY_F32)
		{
			bns::Intersection i2;
			i2.HitShape = this;
			i2.TMinDist = tmax;
			i2.Type = ShapeType::Cube;
			fill_intersections.Add(i2);
		}
	}
}

#pragma endregion

#pragma region CYLINDER SHAPE

bns::CylinderShape::CylinderShape()
{
	// Keep it 1 unit long.
	Minimum = -0.5f;
	Maximum = 0.5f;

	Capped = true;
}

bool bns::CylinderShape::CheckCap(const RayF& ray, F32 t) const 
{
	F32 x = ray.Origin.X + t * ray.Direction.X;
	F32 z = ray.Origin.Z + t * ray.Direction.Z;
	return Pow(x) + Pow(z) <= 1.0f;
}

void bns::CylinderShape::IntersectCaps(const RayF& ray, Intersections& fill_intersections) const
{
	// NOTE: here check if ray.direction.y is close to zero can be used to return and not fill_intersections,
	// but that is already done by calling method, therefore not necessary

	// check for an intesection with the lower end cap by intersecting
	// the ray with the plane at y = minimum
	F32 t = (Minimum - ray.Origin.Y) / ray.Direction.Y;
	if (CheckCap(ray, t))
	{
		Intersection i;
		i.TMinDist = t;
		i.HitShape = this;
		i.Type = ShapeType::Cylinder;
		fill_intersections.Add(i);
	}

	// check for an intesection with the upper end cap by intersecting
	// the ray with the plane at y = maximum
	t = (Maximum - ray.Origin.Y) / ray.Direction.Y;
	if (CheckCap(ray, t))
	{
		Intersection i;
		i.TMinDist = t;
		i.HitShape = this;
		i.Type = ShapeType::Cylinder;
		fill_intersections.Add(i);
	}

}

void bns::CylinderShape::IntersectWithRay(const RayF& ray, Intersections& fill_intersections) const
{
	F32 a = Pow(ray.Direction.X) + Pow(ray.Direction.Z);

	// ray is parallel to y axis, therefore ray is missing the cylinder
	if (Abs(a) < EPSILON) return;

	F32 b = 2.0f * ray.Origin.X * ray.Direction.X + 2.0f * ray.Origin.Z * ray.Direction.Z;
	F32 c = Pow(ray.Origin.X) + Pow(ray.Origin.Z) - 1.0f;

	F32 disc = Pow(b) - 4 * a * c;
	
	// ray does not intersect the cylinder
	if (disc < 0.0f) return;

	disc = Sqrt(disc);
	F32 deno = 2.0f * a;

	F32 t0 = (-b - disc) / deno;
	F32 t1 = (-b + disc) / deno;

	if (t0 > t1) 
	{
		Swap(&t0, &t1);
	}

	F32 y0 = ray.Origin.Y + t0 * ray.Direction.Y;
	if (Minimum < y0 && y0 < Maximum)
	{
		Intersection i0;
		i0.TMinDist = t0;
		i0.HitShape = this;
		i0.Type = ShapeType::Cylinder;
		fill_intersections.Add(i0);
	}

	F32 y1 = ray.Origin.Y + t1 * ray.Direction.Y;
	if (Minimum < y1 && y1 < Maximum)
	{
		Intersection i1;
		i1.TMinDist = t1;
		i1.HitShape = this;
		i1.Type = ShapeType::Cylinder;
		fill_intersections.Add(i1);
	}

	if (this->Capped)
	{
		this->IntersectCaps(ray, fill_intersections);
	}
}

bns::Vec3F bns::CylinderShape::GetLocalNormalAt(bns::Point3F point) const
{
	// sqaure of the distance from y axis
	F32 dist = Pow(point.X) + Pow(point.Z);

	if (dist < 1.0f && point.Y >= Maximum - EPSILON)
	{
		return bns::Vec3F(0.0f, 1.0f, 0.0f);
	}
	else if (dist < 1.0f && point.Y <= Minimum + EPSILON)
	{
		return bns::Vec3F(0.0f, -1.0f, 0.0f);
	}

	return Vec3F(point.X, 0.0f, point.Z);
}

#pragma endregion

#pragma region SHAPE HELPERS

bns::Vec3F bns::GetWorldNormalAt(const bns::BaseShape& shape, bns::Vec3F local_normal)
{
	bns::Vec3F normal = bns::Mat4x4F::Transpose(shape.GetInverseTransform()) * local_normal;
	normal.Normalize();

	if (shape.Parent != nullptr)
	{
		normal = GetWorldNormalAt(*shape.Parent, normal);
	}

	return normal;
}

bns::Point4F bns::WorldPointToLocalPoint(const bns::BaseShape& shape, bns::Point4F world_point)
{
	if (shape.Parent)
	{
		world_point = WorldPointToLocalPoint(*shape.Parent, world_point);
	}

	bns::Point4F result = shape.GetInverseTransform() * world_point;
	return result;
}

#pragma endregion
