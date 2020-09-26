%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

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
function scaleCoeff = LengthUnitToAngstromScale(lengthUnit)
%LengthUnitToAngstromScale() generates the scaling coefficient for
%transforming input length unit to angstrom.
% Input:
%   lengthUnit -- valid input: mm, um, nm, pm
% Output:
%   scaleCoeff -- scaling coefficient, 0 for invalid input

switch lengthUnit
    case 'mm'
        scaleCoeff = 1e7;
    case 'um'
        scaleCoeff = 1e4;
    case 'nm'
        scaleCoeff = 10;
    case 'pm'
        scaleCoeff = 1e-2;
    otherwise
        scaleCoeff = 0;
end

end

