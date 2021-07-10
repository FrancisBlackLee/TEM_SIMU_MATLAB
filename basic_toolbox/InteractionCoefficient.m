function [interCoeff] = InteractionCoefficient(KeV)
%InteractionCofficient.m
%   KeV --electron beam energy (in kilo-electron-volt).

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

lightSpeed = 2.998e8; % in meter
eleCharge = 1.602e-19; % in coulomb
voltage = 1000 * KeV; % in volt
m0 = 9.109e-31; % in kilogram
lambda = HighEnergyWavLen_X(KeV);
interCoeff = 2 * pi/(lambda * KeV) * (m0 * lightSpeed^2 + eleCharge * voltage)/...
    (2 * m0 * lightSpeed^2 + eleCharge * voltage);

end

