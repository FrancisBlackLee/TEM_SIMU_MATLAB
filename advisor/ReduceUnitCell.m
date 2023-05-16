function [periodicUnit, newConvMat] = ReduceUnitCell(unitCell, convMat, maxReduceNum)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2 || nargin > 3
    error('Too few input arguments');
elseif nargin == 2
    maxReduceNum = 10;
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
        for reduceNum = 2 : maxReduceNum
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