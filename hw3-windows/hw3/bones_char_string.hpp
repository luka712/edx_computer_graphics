#pragma once

#ifndef BONES_CHAR_STR_H

#define BONES_CHAR_STR_H

#include "bones_types.hpp"

/// <summary>
/// Finds the integer from string.
/// </summary>
/// <param name="string">String to search.</param>
/// <param name="start_index">Star index in string.</param>
/// <param name="out_end_index">Out parameters of end index, where string search has ended.</param>
/// <param name="out_integer">Result if there is any. 0 if non found, but can also mean that 0 was encountered.</param>
void ReadI32FromString(const char* string, I32 start_index, I32* out_end_index, I32* out_integer);

/// <summary>
/// Finds the integer from string.
/// </summary>
/// <param name="string">String to search.</param>
/// <param name="start_index">Star index in string.</param>
/// <param name="out_integer">Result if there is any. 0 if non found, but can also mean that 0 was encountered.</param>
void ReadI32FromString(const char* string, I32 start_index, I32* out_integer);


/// <summary>
/// Finds the float from string.
/// </summary>
/// <param name="string">String to search.</param>
/// <param name="start_index">Star index in string.</param>
/// <param name="out_end_index">Out parameters of end index, where string search has ended.</param>
/// <param name="out_integer">Result if there is any. 0 if non found, but can also mean that 0 was encountered.</param>
void ReadF32FromString(const char* string, I32 start_index, I32* out_end_index, F32* out_real);

#endif // !BONES_CHAR_STR_H