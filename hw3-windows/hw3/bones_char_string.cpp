#include "bones_char_string.hpp"
#include <string>

void ReadI32FromString(const char* string, I32 start_index, I32* out_end_index, I32* out_integer)
{
	I32 result = 0;

	// increment index until first digit is found
	while (!isdigit(string[start_index]))
	{
		start_index++;
	}

	I32 multiplier = 1;
	if (string[start_index - 1] == '-')
	{
		multiplier = -1;
	}

	while (isdigit(string[start_index]))
	{
		result *= 10;
		result += string[start_index++] - '0';
	}

	*out_end_index = start_index;
	*out_integer = result * multiplier;
}

void ReadI32FromString(const char* string, I32 start_index, I32* out_integer)
{
	// Don't care about end index.
	I32 end_index = 0;
	ReadI32FromString(string, start_index, &end_index, out_integer);
}

void ReadF32FromString(const char* string, I32 start_index, I32* out_end_index, F32* out_real)
{
	F32 result = 0;

	// increment index until first digit is found
	while (!isdigit(string[start_index]))
	{
		start_index++;
	}

	I32 multiplier = 1;
	if (string[start_index - 1] == '-')
	{
		multiplier = -1;
	}

	while (isdigit(string[start_index]))
	{
		result *= 10;
		result += string[start_index++] - '0';
	}

	*out_end_index = start_index;
	*out_real = result * multiplier;
}
