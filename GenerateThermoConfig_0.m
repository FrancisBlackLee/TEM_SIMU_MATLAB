%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outllok.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = GenerateThermoConfig_0(SavePath, StdLatt, MassNum, DebyeTemp, SimuTemp, ConfigNum)
%GenerateThermoConfig_0.m generates a given number of lattice
%configurations under input temperature using Einstein model.
%   SavePath -- save all the generated configuration under a given path, to
%       avoid memory overflow;
%   StdLatt -- standard lattice matrix without displacement, matrix format
%       is [T1, ..., TN; P1, ..., PN; X1, ..., XN; Y1, ..., YN; Z1, ...,
%       ZN], where T denotes atomic type, represented by the atomic
%       numbers, P denotes the elemental proportion, X, Y and Z denote the
%       orthogonal coordinates;
%   MassNum -- mass number array corresponding to the lattice matrix;
%   DebyeTemp -- Debye temperature array corresponding to the lattice
%       matrix;
%   SimuTemp -- simulation temperature;
%   ConfigNum -- number of configurations to be generated;

ThermoDisp = ThermoDisplace_0(MassNum, DebyeTemp, SimuTemp);

for ConfigIdx = 1 : ConfigNum
    % normally distributed random radial displacement:
    rng('shuffle');
    RndRadius = normrnd(zeros(size(ThermoDisp)), ThermoDisp);
    % uniformly distributed random azimuthal angles:
    rng('shuffle');
    RndAzi = pi * rand(size(ThermoDisp));
    % uniformly distributed random polar angles:
    rng('shuffle');
    RndPol = 2 * pi * rand(size(ThermoDisp));

    Xdisp = RndRadius .* sin(RndAzi) .* cos(RndPol);
    Ydisp = RndRadius .* sin(RndAzi) .* sin(RndPol);
    Zdisp = RndRadius .* cos(RndAzi);
    
    DispLatt = StdLatt;
    DispLatt(3, : ) = DispLatt(3, : ) + Xdisp;
    DispLatt(4, : ) = DispLatt(4, : ) + Ydisp;
    DispLatt(5, : ) = DispLatt(5, : ) + Zdisp;
    
    DispLatt = DispLatt';
    filename = strcat(SavePath, 'Config_', num2str(ConfigIdx), '.txt');
    save(filename, 'DispLatt', '-ascii', '-double', '-tabs');
end

end

