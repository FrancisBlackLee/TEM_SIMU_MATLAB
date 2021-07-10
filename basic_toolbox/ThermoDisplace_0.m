function [u] = ThermoDisplace_0(massNum, debyeTemp, temperature)
%ThermoDisplace_0.m calculates the standard deviation of
%thermo-displacement using Einstein model.
%   Input:
%       massNum -- mass number, i.e. for carbon is 12.01;
%       debyeTemp -- Debye temperature (in Kelvin);
%       temperature -- simulation temperature;
%   Output:
%       u -- standard deviation of the thermo-displacement (in
%       Angstrom);

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

u = sqrt(144.38 * temperature ./ (massNum .* debyeTemp.^2));

end

