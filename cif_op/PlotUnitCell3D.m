function PlotUnitCell3D(convMat, fracCoords)
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

% sort fracCoords by the atomic number Z:
[~, atomTypeOrder] = sort(fracCoords(1, :), 'ascend');
fracCoords = fracCoords(:, atomTypeOrder);

atomNum = size(fracCoords, 2);
atomHead = 1;
hold on;
while atomHead <= atomNum
    plotCoords = fracCoords(3 : 5, fracCoords(1, :) == fracCoords(1, atomHead));
    plotCoords = convMat * plotCoords;
    scatter3(plotCoords(1, :), plotCoords(2, :), plotCoords(3, :), 75, 'filled');
    atomHead = atomHead + size(plotCoords, 2);
end

% find all unit cells:
[minFracX, maxFracX] = bounds(fracCoords(3, :));
[minFracY, maxFracY] = bounds(fracCoords(4, :));
[minFracZ, maxFracZ] = bounds(fracCoords(5, :));

lowerBoundA = floor(minFracX);
upperBoundA = ceil(maxFracX);

lowerBoundB = floor(minFracY);
upperBoundB = ceil(maxFracY);

lowerBoundC = floor(minFracZ);
upperBoundC = ceil(maxFracZ);

for aIdx = lowerBoundA : upperBoundA - 1
    for bIdx = lowerBoundB : upperBoundB - 1
        for cIdx = lowerBoundC : upperBoundC - 1
            tmpFracCoords = fracCoords(:, (fracCoords(3, :) >= aIdx) &...
                (fracCoords(3, :) <= aIdx + 1) &...
                (fracCoords(4, :) >= bIdx) &...
                (fracCoords(4, :) <= bIdx + 1) &...
                (fracCoords(5, :) >= cIdx) &...
                (fracCoords(5, :) <= cIdx + 1));
            if ~isempty(tmpFracCoords)
                PlotSingleUnitCell([aIdx, bIdx, cIdx]);
            end
        end
    end
end

hold off;

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');

% nested function:
    function PlotSingleUnitCell(cellIndices)
        initVertexSites = [0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1]';
        adjVertexSites1 = [0, 0, 1; 0, 0, 0; 0, 1, 1; 0, 1, 0; 1, 0, 1; 1, 0, 0; 1, 1, 1; 1, 1, 0]';
        adjVertexSites2 = [0, 1, 0; 0, 1, 1; 0, 0, 0; 0, 0, 1; 1, 1, 0; 1, 1, 1; 1, 0, 0; 1, 0, 1]';
        adjVertexSites3 = [1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1; 0, 0, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1]';
        
        initVertexSites(1, :) = initVertexSites(1, :) + cellIndices(1);
        initVertexSites(2, :) = initVertexSites(2, :) + cellIndices(2);
        initVertexSites(3, :) = initVertexSites(3, :) + cellIndices(3);
        
        adjVertexSites1(1, :) = adjVertexSites1(1, :) + cellIndices(1);
        adjVertexSites1(2, :) = adjVertexSites1(2, :) + cellIndices(2);
        adjVertexSites1(3, :) = adjVertexSites1(3, :) + cellIndices(3);
        
        adjVertexSites2(1, :) = adjVertexSites2(1, :) + cellIndices(1);
        adjVertexSites2(2, :) = adjVertexSites2(2, :) + cellIndices(2);
        adjVertexSites2(3, :) = adjVertexSites2(3, :) + cellIndices(3);
        
        adjVertexSites3(1, :) = adjVertexSites3(1, :) + cellIndices(1);
        adjVertexSites3(2, :) = adjVertexSites3(2, :) + cellIndices(2);
        adjVertexSites3(3, :) = adjVertexSites3(3, :) + cellIndices(3);

        initCellVertexes = convMat * initVertexSites;
        adjCellVertexes1 = convMat * adjVertexSites1;
        adjCellVertexes2 = convMat * adjVertexSites2;
        adjCellVertexes3 = convMat * adjVertexSites3;

        vertexNum = size(initCellVertexes, 2);
        for vertexIdx = 1 : vertexNum
            line([initCellVertexes(1, vertexIdx), adjCellVertexes1(1, vertexIdx)],...
                [initCellVertexes(2, vertexIdx), adjCellVertexes1(2, vertexIdx)],...
                [initCellVertexes(3, vertexIdx), adjCellVertexes1(3, vertexIdx)], 'Color', 'black');

            line([initCellVertexes(1, vertexIdx), adjCellVertexes2(1, vertexIdx)],...
                [initCellVertexes(2, vertexIdx), adjCellVertexes2(2, vertexIdx)],...
                [initCellVertexes(3, vertexIdx), adjCellVertexes2(3, vertexIdx)], 'Color', 'black');

            line([initCellVertexes(1, vertexIdx), adjCellVertexes3(1, vertexIdx)],...
                [initCellVertexes(2, vertexIdx), adjCellVertexes3(2, vertexIdx)],...
                [initCellVertexes(3, vertexIdx), adjCellVertexes3(3, vertexIdx)], 'Color', 'black');
        end
        
    end

end

