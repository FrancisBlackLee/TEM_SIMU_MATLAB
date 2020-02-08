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
% FrozenPhonon_test_1.m -- test the configurations generated in FrozenPhonon_test_0.m
clc;
close all;
clear all;
%% Load configurations:
LoadPath = 'D:\Francis. B. Lee\Practice\Tomography\nvstgt\Si_110_diff_data\300K_configs\';
ConfigNum = 20;

for ConfigIdx = 1 : ConfigNum
    filename = strcat(LoadPath, 'Config_', num2str(ConfigIdx), '.txt');
    TempConfig = load(filename);
    Config{ConfigIdx} = TempConfig;
    figure;
    scatter(TempConfig( : , 4), TempConfig( : , 5));
end