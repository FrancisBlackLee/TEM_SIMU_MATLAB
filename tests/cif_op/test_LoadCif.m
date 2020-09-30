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
% test_LoadCif.m
clc;
clear;
close all;
%% main:
filename_1 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_computed.cif'];

filename_2 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_conventional_standard.cif'];

filename_3 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_primitive.cif'];

filename_4 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_symmetrized.cif'];

filename_5 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'Langbeinite (ortho).cif'];

crysInfo = LoadCif(filename_4);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
atomSiteMat = ExtractAtomSiteFromCrysInfo(crysInfo);
fullAtomSiteMat = AddEquivAtomSites(atomSiteMat);

hkl = [1, 1, 1];
convMat = ConversionMatrix_hkl(cellLengths, cellAngles, hkl);
atomCartCoord = convMat * fullAtomSiteMat(3 : 5, :);

PlotCell2D(convMat, fullAtomSiteMat);
PlotCell3D(convMat, fullAtomSiteMat);

%% extra function
function PlotCell2D(convMat, atomSites)

atomCartCoord = convMat * atomSites(3 : 5, :);
figure;
scatter(atomCartCoord(1, :), atomCartCoord(2, :), 75, 'filled');

initVertexSites = [0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1]';
adjVertexSites1 = [0, 0, 1; 0, 0, 0; 0, 1, 1; 0, 1, 0; 1, 0, 1; 1, 0, 0; 1, 1, 1; 1, 1, 0]';
adjVertexSites2 = [0, 1, 0; 0, 1, 1; 0, 0, 0; 0, 0, 1; 1, 1, 0; 1, 1, 1; 1, 0, 0; 1, 0, 1]';
adjVertexSites3 = [1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1; 0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1]';

initCellVertexes = convMat * initVertexSites;
adjCellVertexes1 = convMat * adjVertexSites1;
adjCellVertexes2 = convMat * adjVertexSites2;
adjCellVertexes3 = convMat * adjVertexSites3;

hold on;
vertexNum = size(initCellVertexes, 2);
for vertexIdx = 1 : vertexNum
    line([initCellVertexes(1, vertexIdx), adjCellVertexes1(1, vertexIdx)],...
        [initCellVertexes(2, vertexIdx), adjCellVertexes1(2, vertexIdx)], 'Color', 'blue');
    
    line([initCellVertexes(1, vertexIdx), adjCellVertexes2(1, vertexIdx)],...
        [initCellVertexes(2, vertexIdx), adjCellVertexes2(2, vertexIdx)], 'Color', 'blue');
    
    line([initCellVertexes(1, vertexIdx), adjCellVertexes3(1, vertexIdx)],...
        [initCellVertexes(2, vertexIdx), adjCellVertexes3(2, vertexIdx)], 'Color', 'blue');
end
hold off;

axis equal;

end

function PlotCell3D(convMat, atomSites)

atomCartCoord = convMat * atomSites(3 : 5, :);
figure;
scatter3(atomCartCoord(1, :), atomCartCoord(2, :), atomCartCoord(3, :), 75, 'filled');

initVertexSites = [0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1]';
adjVertexSites1 = [0, 0, 1; 0, 0, 0; 0, 1, 1; 0, 1, 0; 1, 0, 1; 1, 0, 0; 1, 1, 1; 1, 1, 0]';
adjVertexSites2 = [0, 1, 0; 0, 1, 1; 0, 0, 0; 0, 0, 1; 1, 1, 0; 1, 1, 1; 1, 0, 0; 1, 0, 1]';
adjVertexSites3 = [1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1; 0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1]';

initCellVertexes = convMat * initVertexSites;
adjCellVertexes1 = convMat * adjVertexSites1;
adjCellVertexes2 = convMat * adjVertexSites2;
adjCellVertexes3 = convMat * adjVertexSites3;

hold on;
vertexNum = size(initCellVertexes, 2);
for vertexIdx = 1 : vertexNum
    line([initCellVertexes(1, vertexIdx), adjCellVertexes1(1, vertexIdx)],...
        [initCellVertexes(2, vertexIdx), adjCellVertexes1(2, vertexIdx)],...
        [initCellVertexes(3, vertexIdx), adjCellVertexes1(3, vertexIdx)], 'Color', 'blue');
    
    line([initCellVertexes(1, vertexIdx), adjCellVertexes2(1, vertexIdx)],...
        [initCellVertexes(2, vertexIdx), adjCellVertexes2(2, vertexIdx)],...
        [initCellVertexes(3, vertexIdx), adjCellVertexes2(3, vertexIdx)], 'Color', 'blue');
    
    line([initCellVertexes(1, vertexIdx), adjCellVertexes3(1, vertexIdx)],...
        [initCellVertexes(2, vertexIdx), adjCellVertexes3(2, vertexIdx)],...
        [initCellVertexes(3, vertexIdx), adjCellVertexes3(3, vertexIdx)], 'Color', 'blue');
end
hold off;

axis equal;

end