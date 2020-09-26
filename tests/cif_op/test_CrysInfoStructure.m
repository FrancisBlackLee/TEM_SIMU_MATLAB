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
% test_CrysInfoStructure.m
clc;
clear;
close all;
%% main:
for i = 1 : 5
    strCellArray = {num2str(i), num2str(i + 1)};
    crysInfo.valuedProperty{i, 1} = strCellArray{1};
    crysInfo.valuedProperty{i, 2} = strCellArray{2};
end

for i = 1 : 5
    strCellArray = num2str(i);
    crysInfo.loopProperty{i, 1} = strCellArray;
end

data = {'1', '2', '3', '4', '5'};

for i = 1 : 5
    crysInfo.loopProperty{i, 2} = data{i};
end