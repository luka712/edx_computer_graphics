#include "bones_lights.hpp"

bns::BaseLight::BaseLight(LightType type)
	:Type(type)
{

}

bns::DirectionalLight::DirectionalLight()
	:BaseLight(LightType::Directional)
{
}

bns::PointLight::PointLight()
:BaseLight(LightType::Point)
{
	Attenuation.Constant = 1.0f;
	Attenuation.Linear = 0.0f;
	Attenuation.Quadratic = 0.0f;
}
