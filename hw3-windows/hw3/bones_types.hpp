#pragma once
#ifndef BONES_TYPES_H

#define BONES_TYPES_H

/*
	All of the core types for every "bones" library component, structure, class, function.
	To be used by all of the "bones" related libraries.
*/

#include <cstdint>
#include <limits>


namespace bns
{
	typedef float F32;
	typedef uint8_t U8;
	typedef uint32_t U32;
	typedef int32_t I32;

	const static F32 MAX_F32 = std::numeric_limits<F32>::max();
	const static F32 INFINITY_F32 = std::numeric_limits<F32>::infinity();
	const static F32 NEGATIVE_INFINITY_F32 = std::numeric_limits<F32>::infinity() * -1.0f;
}

#endif // !BONES_TYPES_H

