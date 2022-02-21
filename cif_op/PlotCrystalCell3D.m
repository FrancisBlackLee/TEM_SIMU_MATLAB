function PlotCrystalCell3D(convMat, fracCoords)
%PlotCrystalCell2D.m plots the atoms in 3D cartesian coordinates, given the
%conversion matrix and fractional coordinates.
%   convMat -- conversion matrix;
%   fracCoords -- fractional coordinates of the unit cell, syntax: [T; P;
%       fracX; fracY; fracZ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2022  Francis Black Lee (Li Xian)

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

fracCoords = AddEquivAtomSites(fracCoords);
atomCartCoord = convMat * fracCoords(3 : 5, :);
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
xlabel('x');
ylabel('y');
zlabel('z');

end

