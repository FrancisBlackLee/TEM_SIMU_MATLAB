/*
	status.h is the open source part of vtemlab v0.0 engine,
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

#ifndef STATUS_H
#define STATUS_H

enum
{
	VTEMLAB_SUCCESS  = 0,  // success
	VTEMLAB_FAILURE  = -1, // failure
	VTEMLAB_ENULLPTR = 1,  // null pointer
	VTEMLAB_ENOMEM   = 2,  // malloc failure, no memory
	VTEMLAB_ELOSS    = 3,  // loss of accuracy
	VTEMLAB_EFILEW   = 4,  // file writing error
	VTEMLAB_EFILER   = 5,  // file reading error
	VTEMLAB_ENOENT   = 6,  // directory does not exist
	VTEMLAB_EINVAL   = 7   // invalid argument 
};

// print error message for error object:
void PrintErrorMsg(const char* errorObj, const int errorCode);

// Update process:
void UpdateProcess(const char* processName, int doneNum, int totalNum);

#endif // !STATUS_H
