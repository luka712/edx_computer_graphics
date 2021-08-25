#pragma once
#include "bones_types.hpp"

#ifndef BONES_FILES_H

#define BONES_FILES_H

namespace bns
{
	struct FileContents
	{
		char* Contents;
		U32 Length;

		~FileContents();
	};

	FileContents* ReadAndCloseFile(const char* filename);
	void FreeFileContents(FileContents* contents);

}

#endif // !BONES_FILES_H

