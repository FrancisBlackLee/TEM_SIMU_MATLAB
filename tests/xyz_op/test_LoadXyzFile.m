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
% test_LoadXyzFile.m
clc;
clear;
close all;
%% main:
filename = 'Tests\Nanocluster.xyz';
filePtr = fopen(filename, 'r');

atomNum = fscanf(filePtr, '%d', 1);
atomTypeList = zeros(1, atomNum);
atomXyzList = zeros(3, atomNum);
for atomIdx = 1 : atomNum
    tmpTypeStr = fscanf(filePtr, '%s', 1);
    atomTypeList(atomIdx) = AtomTypeStrToIdx(tmpTypeStr);
    atomXyzList(1, atomIdx) = fscanf(filePtr, '%f', 1);
    atomXyzList(2, atomIdx) = fscanf(filePtr, '%f', 1);
    atomXyzList(3, atomIdx) = fscanf(filePtr, '%f', 1);
end

fclose(filePtr);

figure;
subplot(2, 2, 1);
scatter3(atomXyzList(1, :), atomXyzList(2, :), atomXyzList(3, :), '.');
title('3D');
axis square;

subplot(2, 2, 2);
scatter(atomXyzList(1, :), atomXyzList(2, :), '.');
title('x-y');
axis square;

subplot(2, 2, 3);
scatter(atomXyzList(1, :), atomXyzList(3, :), '.');
title('x-z');
axis square;

subplot(2, 2, 4);
scatter(atomXyzList(2, :), atomXyzList(3, :), '.');
title('y-z');
axis square;

%% Test:
[testTypeList, testXyzList] = ReadCrystalMakerXyz(filename);
if isequal(testTypeList, atomTypeList) && isequal(testXyzList, atomXyzList)
    msgbox('PASSED');
else
    msgbox('FAILED');
end