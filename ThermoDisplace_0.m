%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

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
function [MeanDisplace] = ThermoDisplace_0(MassNum, DebyeTemp, Temp)
%ThermoDisplace_0.m calculates the standard deviation of
%thermo-displacement using Einstein model.
%   Input:
%       MassNum -- mass number, i.e. for carbon is 12.01;
%       DebyeTemp -- Debye temperature (in Kelvin);
%       Temp -- simulation temperature;
%   Output:
%       MeanDisplace -- standard deviation of the thermo-displacement (in
%       Angstrom);

MeanDisplace = sqrt(144.38 * Temp ./ (MassNum .* DebyeTemp.^2));

end

