/*
	mode.h is the open source part of vtemlab v0.0 engine,
	providing operation modes for vtemlab.

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

#ifndef MODE_H
#define MODE_H

enum
{
	FFTSHIFTED_VEC      = 101,
	UNFFTSHIFTED_VEC    = 102,
	ENCRYPTED_PROJPOT   = 103,
	DECRYPTED_PROJPOT   = 104,
	CBED_UNLOG          = 105,
	CBED_LOG            = 106,
	FRAC_COORD          = 107,
	CART_COORD          = 108,
	COMPLETE_BOUNDARY   = 109,
	INCOMPLETE_BOUNDARY = 200
};

#endif // !MODE_H
