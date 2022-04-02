/*
	status.cpp is the open source part of vtemlab v0.0 engine,
	providing status info/operations for vtemlab.

	Copyright (C) 2022  Francis Black Lee (Chinese name: ¿ÓœÕ, Li Xian)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.

	Email: warner323@outlook.com
*/

#include <stdio.h>

#include "status.h"

const char* ErrorMessage(const int errorCode)
{
	switch (errorCode)
	{
	case VTEMLAB_SUCCESS:
		return "Succes";
		break;
	case VTEMLAB_FAILURE:
		return "Failure";
		break;
	case VTEMLAB_ENULLPTR:
		return "NULL pointer";
		break;
	case VTEMLAB_ENOMEM:
		return "Memory allocation failed";
		break;
	case VTEMLAB_ELOSS:
		return "Loss of accuracy";
		break;
	case VTEMLAB_EFILER:
		return "File reading failure";
		break;
	case VTEMLAB_EFILEW:
		return "File writing failure";
		break;
	case VTEMLAB_ENOENT:
		return "Directory or file does not exist";
		break;
	case VTEMLAB_EINVAL:
		return "Invalid argument";
		break;
	default:
		return "Undefined error";
		break;
	}
}


// print error message for error object:
void PrintErrorMsg(const char* errorObj, const int errorCode)
{
	const char* errorMsg = ErrorMessage(errorCode);
	printf("%s: %s\n\n", errorObj, errorMsg);
}


// Update process:
void UpdateProcess(const char* processName, int doneNum, int totalNum)
{
	printf("\r");
	printf("%s: %d / %d completed...", processName, doneNum, totalNum);
}