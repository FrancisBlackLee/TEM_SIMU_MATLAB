function [periodicUnit, newConvMat] = ReduceUnitCell(unitCell, convMat, maxIter)
%ReduceUnitCell reduces the input unit cell to a smaller unit cell and
%adjusts the conversion matrix as well.
% Input:
%   unitCell -- unit cell to be reduced;
%   convMat -- conversion matrix of the input unit cell;
%   maxIter -- max iteration along each direction.
% Output:
%   periodicUnit -- reduced unit cell, which may be a smaller periodic
%       unit;
%   newConvMat -- new conversion matrix of the reduced unit cell.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

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

if nargin < 2 || nargin > 3
    error('Too few input arguments');
elseif nargin == 2
    maxIter = 10;
end

fullUnitCell = AddEquivAtomSites(unitCell);
[tmpUnit, nx] = ReduceAlongDim(fullUnitCell, 3);
[tmpUnit, ny] = ReduceAlongDim(tmpUnit, 4);
[tmpUnit, nz] = ReduceAlongDim(tmpUnit, 5);

periodicUnit = tmpUnit;
periodicUnit(3, :) = nx * periodicUnit(3, :);
periodicUnit(4, :) = ny * periodicUnit(4, :);
periodicUnit(5, :) = nz * periodicUnit(5, :);

newConvMat = convMat;
newConvMat(:, 1) = newConvMat(:, 1) / nx;
newConvMat(:, 2) = newConvMat(:, 2) / ny;
newConvMat(:, 3) = newConvMat(:, 3) / nz;

% nested function
    function [outUnit, n] = ReduceAlongDim(inUnit, dim)
        outUnit = inUnit;
        n = 1;
        for reduceNum = 2 : maxIter
            glide = 1 / reduceNum;
            glideUnit = inUnit(:, inUnit(dim, :) <= glide);
            if ~isempty(glideUnit)
                isPeriodic = true;
                for glideIdx = 1 : reduceNum - 1
                    tmpGlidedUnit = glideUnit;
                    tmpGlidedUnit(dim, :) = tmpGlidedUnit(dim, :) + glideIdx * glide;
                    lia = ismembertol(tmpGlidedUnit', inUnit', 'ByRows', true);
        
                    if sum(lia, 'all') ~= size(glideUnit, 2)
                        isPeriodic = false;
                        break;
                    end
                end
        
                if isPeriodic
                    outUnit = glideUnit;
                    n = reduceNum;
                end
            end
        end
    end


end