function [atomCoordMat] = CreateNanoCluster_uvw(atomSiteMat, cellLengths,...
    cellAngles, uvw, radius)
%CreateNanoCluster_uvw() rotates the unit cell to the given orientation, 
%duplicate the unit cell and rehapes the lattice to nano-cluster.
% Input:
%   atomSiteMat -- atomic site matrix;
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases a and c);
%       2 for cell angle beta (between bases b and c) and
%       3 for cell angle gamma (between bases a and b);
%   uvw -- Orientation indices;
%   radius -- radius of the nano-cluster;
% Output:
%   atomCoordMat -- atomic coordinate matrix (cartesian);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

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

if any(mod(uvw, 1))
    errID = 'myComponent:inputError';
    msgtext = 'Orientation indices must be integers.';
    ME = MException(errID, msgtext);
    throw(ME);
elseif isempty(atomSiteMat)
    errID = 'myComponent:inputError';
    msgtext = 'Atomic site matrix cannot be empty.';
    ME = MException(errID, msgtext);
    throw(ME);
elseif ~all(cellLengths)
    errID = 'myComponent:inputError';
    msgtext = 'Cell lengths cannot contain zero.';
    ME = MException(errID, msgtext);
    throw(ME);
elseif ~all(cellAngles)
    errID = 'myComponent:inputError';
    msgtext = 'Cell angles cannot contain zero.';
    ME = MException(errID, msgtext);
    throw(ME);
elseif ~any(uvw)
    errID = 'myComponent:inputError';
    msgtext = 'Orientation indices cannot be all zeros.';
    ME = MException(errID, msgtext);
    throw(ME);
elseif ~(radius > 0)
    errID = 'myComponent:inputError';
    msgtext = 'Radius must have positive value.';
    ME = MException(errID, msgtext);
    throw(ME);
else
    convMat = ConversionMatrix_uvw(cellLengths, cellAngles, uvw);
    vecA = convMat(:, 1);
    vecB = convMat(:, 2);
    vecC = convMat(:, 3);

    areaAB = norm(cross(vecA, vecB));
    areaBC = norm(cross(vecB, vecC));
    areaCA = norm(cross(vecC, vecA));
    volume = abs(dot(vecA, cross(vecB, vecC)));
    diameter = 2 * radius;
    cellNumA = ceil(diameter / (volume / areaBC));
    cellNumB = ceil(diameter / (volume / areaCA));
    cellNumC = ceil(diameter / (volume / areaAB));
    
    repStart = 1;
    initAtomNum = size(atomSiteMat, 2);
    expanNum = cellNumA * cellNumB * cellNumC;
    atomCoordMat = repmat(atomSiteMat, 1, expanNum);
    wbHandle = waitbar(0, 'Processing...');
    totalTask = cellNumA * cellNumB * cellNumC;
    finishedTask = 0;
    for aIdx = 0 : cellNumA - 1
        for bIdx = 0 : cellNumB - 1
            for cIdx = 0 : cellNumC - 1
                atomCoordMat(3, repStart : repStart + initAtomNum - 1) =...
                    atomCoordMat(3, repStart : repStart + initAtomNum - 1) + aIdx;
                atomCoordMat(4, repStart : repStart + initAtomNum - 1) =...
                    atomCoordMat(4, repStart : repStart + initAtomNum - 1) + bIdx;
                atomCoordMat(5, repStart : repStart + initAtomNum - 1) =...
                    atomCoordMat(5, repStart : repStart + initAtomNum - 1) + cIdx;
                
                repStart = repStart + initAtomNum;
                finishedTask = finishedTask + 1;
                waitbar(finishedTask / totalTask, wbHandle, 'Processing...');
            end
        end
    end
                    
    tolerance = 1e-8;
    atomCoordMat = atomCoordMat';
    atomCoordMat = uniquetol(atomCoordMat, tolerance, 'ByRows', true);
    atomCoordMat = atomCoordMat';
    
    atomCoordMat(3 : 5, :) = convMat * atomCoordMat(3 : 5, :);
    center = (cellNumA * convMat(:, 1) + cellNumB * convMat(:, 2) + cellNumC * convMat(:, 3)) / 2;
    atomCoordMat(3 : 5, :) = atomCoordMat(3 : 5, :) - center;
    distToOrigin = sqrt(sum(atomCoordMat(3 : 5, :).^2, 1));
    atomCoordMat(:, distToOrigin > radius) = [];
    
    close(wbHandle);
end

end

