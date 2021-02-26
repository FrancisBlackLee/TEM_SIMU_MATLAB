%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019- 2021  Francis Black Lee and Li Xian

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
function [InterCoeff] = InteractionCoefficient(KeV)
%InteractionCofficient.m
%   KeV --electron beam energy (in kilo-electron-volt).

LightSpeed = 2.998e8; % in meter
EleCharge = 1.602e-19; % in coulomb
V = 1000 * KeV; % in volt
m0 = 9.109e-31; % in kilogram
lambda = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);
InterCoeff = 2 * pi/(lambda * KeV) * (m0 * LightSpeed^2 + EleCharge * V)/(2 * m0 * LightSpeed^2 + EleCharge * V);

end

