function [b1, b2, b3] = DirectToReciprocal(a1, a2, a3)
%DirectToReciprocal calculates the reciprocal lattice vectors using the
%direct lattice vectors.
%   a1, a2, a3 -- direct lattice vectors;
%   b1, b2, b3 -- reciprocal lattice vectors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vc = dot(a1, cross(a2, a3));
b1 = 2 * pi * cross(a2, a3) / vc;
b2 = 2 * pi * cross(a3, a1) / vc;
b3 = 2 * pi * cross(a1, a2) / vc;

end