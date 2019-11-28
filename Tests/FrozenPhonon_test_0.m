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
% FrozenPhonon_test_0.m -- generate several configurations for test:
clc;
close all;
clear all;
%% Si <1 1 0>
LattConst = [3.8396, 5.4300, 3.8396]; % [a b c]
CellNum = [9, 6, 15];
PrimeLatt = [14, 14, 14, 14;...
           1, 1, 1, 1;
           0, 0.5, 0, 0.5;...
           0, 0.75, 0.25, 0.5;...
           0, 0, 0.5, 0.5];
PrimeMassNum = [28.086, 28.086, 28.086, 28.086];
PrimeDebyeTemp = [692, 692, 692, 692];
ExpLatt = SquareLattExpanX(PrimeLatt, LattConst, CellNum, 1e-5, PrimeMassNum, PrimeDebyeTemp);
StdLatt = ExpLatt(1 : 5, : );
MassNum = ExpLatt(6, : );
DebyeTemp = ExpLatt(7, : );
SimuTemp = 300;
ConfigNum = 20;
SavePath = 'D:\Francis. B. Lee\Practice\Tomography\nvstgt\Si_110_diff_data\300K_configs\';

RefStdLatt = StdLatt';
save(strcat(SavePath, 'StdLatt.txt'), 'RefStdLatt', '-ascii', '-double', '-tabs');
GenerateThermoConfig_0(SavePath, StdLatt, MassNum, DebyeTemp, SimuTemp, ConfigNum);