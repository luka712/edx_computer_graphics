#include "bones_files.hpp"
#include <cstdio>
#include <malloc.h>

bns::FileContents::~FileContents()
{
	free(Contents);
}

bns::FileContents* bns::ReadAndCloseFile(const char* filename)
{
	bns::FileContents* result = new bns::FileContents();

	FILE* fp = 0;

	fopen_s(&fp, filename, "r");

	// if there if pointer to file
	if (fp != 0)
	{
		// find end and keep that as length
		fseek(fp, 0, SEEK_END);
		result->Length = ftell(fp);
		fseek(fp, 0, SEEK_SET);

		result->Contents = (char*)malloc(result->Length);
		if (result->Contents != 0)
		{
			fread(result->Contents, 1, result->Length, fp);
		}

		fclose(fp);
	}

	return result;
}

void bns::FreeFileContents(FileContents* contents)
{
	delete contents;
}
