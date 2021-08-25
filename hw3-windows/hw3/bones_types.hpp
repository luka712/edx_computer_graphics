#pragma once
#ifndef BONES_TYPES_H

#define BONES_TYPES_H

/*
	All of the core types for every "bones" library component, structure, class, function.
	To be used by all of the "bones" related libraries.
*/

#include <cstdint>
#include <limits>


// TODO: move to namespace

typedef float F32;
typedef uint32_t U32;
typedef int32_t I32;

const static F32 MAX_F32 = std::numeric_limits<F32>::max();

#endif // !BONES_TYPES_H

