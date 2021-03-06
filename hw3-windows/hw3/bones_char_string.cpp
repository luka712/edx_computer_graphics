#include "bones_char_string.hpp"
#include <string>

bool bns::ReadI32FromString(const char* string, bns::I32 start_index, bns::I32* out_end_index, bns::I32* out_integer)
{
	bns::I32 result = 0;

	if (string[start_index] == '\0')
	{
		return false;
	}

	// increment index until first digit is found
	while (!isdigit(string[start_index]))
	{
		start_index++;
	}

	bns::I32 multiplier = 1;
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

	return true;
}

bool bns::ReadI32FromString(const char* string, bns::I32 start_index, bns::I32* out_integer)
{
	// Don't care about end index.
	bns::I32 end_index = 0;
	return ReadI32FromString(string, start_index, &end_index, out_integer);
}

void bns::ReadF32FromString(const char* string, bns::I32 start_index, bns::I32* out_end_index, bns::F32* out_real)
{
	bns::F32 result = 0;

	// increment index until first digit is found, or '.' is found with digit after dot
	while (
		!(
			isdigit(string[start_index]) ||
			(string[start_index] == '.' && isdigit(string[start_index+1]))
		)
	)
	{
		start_index++;
	}

	bns::I32 multiplier = 1;
	if (string[start_index - 1] == '-')
	{
		multiplier = -1;
	}

	
	bool is_decimal = false;
	bns::F32 denominator = 0.0f;
	while (isdigit(string[start_index]) || string[start_index] == '.')
	{
		if (string[start_index] == '.')
		{
			is_decimal = true;
		}
		else 
		{
			result *= 10.0f;
			result += string[start_index] - '0';
			if (is_decimal)
			{
				if (denominator == 0.0f)
				{
					denominator = 10.0f;
				}
				else
				{
					denominator *= 10.0f;
				}
			}
		}

		start_index++;
	}

	if (is_decimal)
	{
		result /= denominator;
	}


	*out_end_index = start_index;
	*out_real = result * multiplier;
}
