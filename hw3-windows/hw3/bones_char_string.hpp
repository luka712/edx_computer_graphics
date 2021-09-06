#pragma once

#ifndef BONES_CHAR_STR_H

#define BONES_CHAR_STR_H

#include "bones_types.hpp"

namespace bns
{

	/// <summary>
	/// Finds the integer from string.
	/// While there is something to read returns true.
	/// </summary>
	/// <param name="string">String to search.</param>
	/// <param name="start_index">Star index in string.</param>
	/// <param name="out_end_index">Out parameters of end index, where string search has ended.</param>
	/// <param name="out_integer">Result if there is any. 0 if non found, but can also mean that 0 was encountered.</param>
	bool ReadI32FromString(const char* string, bns::I32 start_index, bns::I32* out_end_index, bns::I32* out_integer);

	/// <summary>
	/// Finds the integer from string.
	/// While there is something to read returns true.
	/// </summary>
	/// <param name="string">String to search.</param>
	/// <param name="start_index">Star index in string.</param>
	/// <param name="out_integer">Result if there is any. 0 if non found, but can also mean that 0 was encountered.</param>
	bool ReadI32FromString(const char* string, bns::I32 start_index, bns::I32* out_integer);


	/// <summary>
	/// Finds the float from string.
	/// </summary>
	/// <param name="string">String to search.</param>
	/// <param name="start_index">Star index in string.</param>
	/// <param name="out_end_index">Out parameters of end index, where string search has ended.</param>
	/// <param name="out_integer">Result if there is any. 0 if non found, but can also mean that 0 was encountered.</param>
	void ReadF32FromString(const char* string, bns::I32 start_index, bns::I32* out_end_index, bns::F32* out_real);

}

#endif // !BONES_CHAR_STR_H
