#include "bones_raytracer.hpp"

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

	bns::RayF ray(eye, dir);
	return ray;
}

bns::ColorF bns::FindColor(Intersection hit)
{

	// TODO: implement
	if (hit.MinDist != MAX_F32)
	{
		return hit.HitShape->Material.Color;
	}
	return bns::ColorF(0, 0, 0);
}

bns::Intersection bns::GetIntersection(const bns::RayF& ray, bns::BaseShape** shapes, U32 Count)
{
	F32 min_dist = MAX_F32;
	const BaseShape* ptr_hit_shape = nullptr;
	for (U32 i = 0; i < Count; i++)
	{
		const BaseShape* shape = shapes[i];
		F32Pair t = shape->IntersectionDistance(ray);

		// if both t1 and t2 are being less then 0, discard as this ray doesn't intersect anything
		if (t.t1 < 0 || t.t2 < 0)
		{
			continue;
		}
		
		// check if any pair is lower then distance
		if (t.t1 < min_dist || t.t2 < min_dist)
		{
			// if so select lower one 
			min_dist = t.t1 < t.t2 ? t.t1 : t.t2;
			ptr_hit_shape = shapes[i];
		}
	}

	Intersection result;
	result.MinDist = min_dist;
	result.HitShape = ptr_hit_shape;
	result.Type = ShapeType::Triangle;

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

void bns::RayTrace(const Camera& cam, bns::BaseShape** shapes, U32 count_of_shapes, bns::ColorF** colors_to_fill)
{
	for (size_t j = 0; j < cam.ScreenHeight; j++)
	{
		for (size_t i = 0; i < cam.ScreenWidth; i++)
		{
			bns::RayF ray = bns::RayThroughPixel(cam, i, j);
			bns::Intersection hit = GetIntersection(ray, shapes, count_of_shapes);
			colors_to_fill[j][i] = bns::FindColor(hit);
		}
	}
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
}


bns::TriangleShape::TriangleShape()
{
	Triangle = bns::TriangleF();
}

bns::F32Pair bns::TriangleShape::IntersectionDistance(const bns::RayF& ray) const
{
	F32 result = ray.IntersectionDistanceWithTriangle(this->Triangle);
	return
	{
		result, result
	};
}

bns::F32Pair bns::SphereShape::IntersectionDistance(const bns::RayF& ray) const
{
	F32 t1, t2;
	ray.IntersectionDistanceWithSphere(this->Sphere, &t1, &t2);

	return  { t1, t2 };
}

bns::SphereShape::SphereShape()
{
	Sphere = bns::SphereF();
}
