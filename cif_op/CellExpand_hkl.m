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
function [atomCoordMat] = CellExpand_hkl(atomSiteMat, cellLengths,...
    cellAngles, hkl, sideLength)
%CellExpand_hkl() rotates the unit cell to the given orientation, duplicate
%the unit cell and rehapes the cell.
% Input:
%   atomSiteMat -- atomic site matrix;
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases a and c);
%       2 for cell angle beta (between bases b and c) and
%       3 for cell angle gamma (between bases a and b);
%   hkl -- Miller indices;
%   sideLength -- the size of the square sample you want;
% Output:
%   atomCoordMat -- atomic coordinate matrix (cartesian);

if (~any(mod(hkl, 1))) && (~isempty(atomSiteMat)) && all(cellLengths) && all(cellAngles) && any(hkl)
    convMat = ConversionMatrix_hkl(cellLengths, cellAngles, hkl);
    vecZ = [0; 0; 1];
    cosVecAZ = max(min(dot(convMat(:, 1), vecZ) /...
        (norm(convMat(:, 1)) * norm(vecZ)), 1), -1);
    cosVecBZ = max(min(dot(convMat(:, 2), vecZ) /...
        (norm(convMat(:, 2)) * norm(vecZ)), 1), -1);
    cosVecCZ = max(min(dot(convMat(:, 3), vecZ) /...
        (norm(convMat(:, 3)) * norm(vecZ)), 1), -1);
    
    [maxCos, maxCosIdx] = max(abs([cosVecAZ, cosVecBZ, cosVecCZ]));
    repBasisIndices = 1 : 3;
    repBasisIndices(maxCosIdx) = [];
    
    projBasis_1 = [convMat(1, repBasisIndices(1)); convMat(2, repBasisIndices(1))];
    projBasis_2 = [convMat(1, repBasisIndices(2)); convMat(2, repBasisIndices(2))];
    
    cosProjBasisAngle = max(min(dot(projBasis_1, projBasis_2) /...
        (norm(projBasis_1) * norm(projBasis_2)), 1), -1);
    sinProjBasisAngle = sqrt(1 - cosProjBasisAngle^2);
    if sinProjBasisAngle > 0
        % 1.5 is magic number
        projBasis_1_repNum = ceil(sqrt(2) * 1.5 * sideLength / (norm(projBasis_1) * sinProjBasisAngle));
        projBasis_2_repNum = ceil(sqrt(2) * 1.5 * sideLength / (norm(projBasis_2) * sinProjBasisAngle));
        
        newHkl = hkl;
        newHkl(repBasisIndices(1)) = projBasis_1_repNum;
        newHkl(repBasisIndices(2)) = projBasis_2_repNum;
        
        repeatNum = 1;
        for i = 1 : 3
            if newHkl(i) ~= 0
                repeatNum = repeatNum * newHkl(i);
            end
        end
        repAtomSiteMat = repmat(atomSiteMat, 1, repeatNum);
        repSigns = sign(newHkl);

        initAtomNum = size(atomSiteMat, 2);
        repStart = 1;
        for hIdx = 0 : (abs(newHkl(1)) > 1) * (abs(newHkl(1)) - 1)
            for kIdx = 0 : (abs(newHkl(2)) > 1) * (abs(newHkl(2)) - 1)
                for lIdx = 0 : (abs(newHkl(3)) > 1) * (abs(newHkl(3)) - 1)
                    repAtomSiteMat(3, repStart : repStart + initAtomNum - 1) =...
                        repAtomSiteMat(3, repStart : repStart + initAtomNum - 1) + repSigns(1) * hIdx;
                    repAtomSiteMat(4, repStart : repStart + initAtomNum - 1) =...
                        repAtomSiteMat(4, repStart : repStart + initAtomNum - 1) + repSigns(2) * kIdx;
                    repAtomSiteMat(5, repStart : repStart + initAtomNum - 1) =...
                        repAtomSiteMat(5, repStart : repStart + initAtomNum - 1) + repSigns(3) * lIdx;
                    
                    repStart = repStart + initAtomNum;
                end
            end
        end
        
        tolerance = 1e-8;
        
        repAtomSiteMat = repAtomSiteMat';
        repAtomSiteMat = uniquetol(repAtomSiteMat, tolerance, 'ByRows', true);
        repAtomSiteMat = repAtomSiteMat';
        
        atomCoordMat = repAtomSiteMat;
        atomCoordMat(3 : 5, :) = convMat * atomCoordMat(3 : 5, :);
        center = mean(atomCoordMat(3 : 4, :), 2);
        atomCoordMat(3 : 4, :) = atomCoordMat(3 : 4, :) - center;
        atomCoordMat(5, :) = atomCoordMat(5, :) - min(atomCoordMat(5, :));
        
        projVecZ = hkl(1) * convMat(:, 1) + hkl(2) * convMat(:, 2) + hkl(3) * convMat(:, 3);
        projVecZNorm = norm(projVecZ);
        atomCoordMat(5, :) = mod(atomCoordMat(5, :), projVecZNorm);
        
        atomCoordMat = atomCoordMat';
        atomCoordMat = uniquetol(atomCoordMat, tolerance, 'ByRows', true);
        atomCoordMat = atomCoordMat';
    else
        msgbox('Invalid input!');
        atomCoordMat = 0;
    end
else
    msgbox('Invalid input!');
    atomCoordMat = 0;
end


end

