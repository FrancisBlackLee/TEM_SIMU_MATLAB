function scaleCoeff = AngleUnitToRadianScale(angleUnit)
%AngleUnitToRadianScale() generates the scaling coefficient for converting
%the angle unit to radian.
% Input:
%   angleUnit -- valid value: 'mrad', 'rad', 'degree';
% Output:
%   scaleCoeff -- scaling coeffcient, 0 for invalid input;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

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

switch angleUnit
    case 'mrad'
        scaleCoeff = 1.0e-3;
    case 'rad'
        scaleCoeff = 1.0;
    case 'degree'
        scaleCoeff = pi / 180;
    otherwise
        scaleCoeff = 0;
end

end

