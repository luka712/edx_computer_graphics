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
}
