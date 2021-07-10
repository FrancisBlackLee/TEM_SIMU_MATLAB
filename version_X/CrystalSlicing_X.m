function [slice, sliceDist, extraSlice] = CrystalSlicing_X(stdLatt, ...
    thermoLatt, maxSliceSpacing, zMax, lattMode, YN, plotColor)
%CrystalSlicing_X.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   StdLatt -- standard lattice matrix (without thermo-displacements), uses
%       the most recently accepted crystal info matrix format: [T; P; X; Y;
%       Z], T: atomic types (atomic numbers); P: elemental proportion; X,
%       Y and Z: orthogonal atomic coordinates.
%   ThermoLatt -- lattice with thermo-displacements, adopt the same format
%       of StdLatt;
%   DistError -- the largest error distance to judge whether atoms of
%       different heights be rearranged to one slice;
%   YN -- whether to show each slice: 1 --yes, 0 --no.
%   PlotColor -- colors for each type of atom;
%   LattMode -- lattice mode:
%       LattMode == 0 : cubic specimen;
%       LattMode == 1: nanoparticle;
%   NOTE: X denotes an experimental version!
%   Current target: lattice with thermo-displacements.

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

atomNum = size(stdLatt, 2);
if nargin == 6
    plotColor = ones(1, atomNum);
end

[zStd, zOrder] = sort(stdLatt(5, : ), 'descend');
stdLatt = stdLatt( : , zOrder);
thermoLatt = thermoLatt( : , zOrder);
plotColor = plotColor(zOrder);

sliceNum = 1;
slice{sliceNum} = thermoLatt( : , 1);
color{sliceNum} = plotColor(1);
tempZ = zStd(1);
% disp(TempZ);
for zIdx = 2 : length(zStd)
    if abs(zStd(zIdx) - tempZ) >= maxSliceSpacing
        sliceDist(sliceNum) = abs(zStd(zIdx) - tempZ);
        sliceNum = sliceNum + 1;
        tempZ = zStd(zIdx);
%         disp(TempZ);
        thermoLatt(5, zIdx) = tempZ;
        slice{sliceNum} = thermoLatt( : , zIdx);
        color{sliceNum} = plotColor(zIdx);
    else
        thermoLatt(5, zIdx) = tempZ;
        slice{sliceNum} = [slice{sliceNum}, thermoLatt( : , zIdx)];
        color{sliceNum} = [color{sliceNum}, plotColor(zIdx)];
    end
end
% determine whether the last slice should be preserved, deleted or counted
% as the first slice:
if lattMode == 1
    sliceDist(sliceNum) = zMax - sum(sliceDist);
    extraSlice = 0;
    xmin = min(stdLatt(3, : )) - 1;
    xmax = max(stdLatt(3, : )) + 1;
    ymin = min(stdLatt(4, : )) - 1;
    ymax = max(stdLatt(4, : )) + 1;
    if YN == 1
        for SliceIdx = 1 : sliceNum
            ColorList = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'];
            figure;
            hold on;
            for AtomIdx = 1 : size(slice{SliceIdx}, 2)
                scatter(slice{SliceIdx}(3, AtomIdx), slice{SliceIdx}(4, AtomIdx), 'o', ColorList(color{SliceIdx}(AtomIdx)));
            end
            axis([xmin xmax ymin ymax]);
            axis equal;
            title(['slice ', num2str(SliceIdx)]);
        end
    end
else
    if zMax - sum(sliceDist) >= maxSliceSpacing
        sliceDist(sliceNum) = zMax - sum(sliceDist);
        extraSlice = 0;
    else
        % rearrange the first slice and last slice by ascending order of x,
        % then y, so that whether the last slice is identical to the first
        % slice could be determined
        IdMinDist = 1e-2;
        IdCrit = 0;
        for TopIdx = 1 : size(slice{1}, 2)
            for BotIdx = 1 : size(slice{sliceNum}, 2)
                TempDist = sqrt(sum((slice{1}(3:4, TopIdx) - slice{sliceNum}(3:4, BotIdx)).^2));
                if (TempDist <= IdMinDist) && (slice{1}(1, TopIdx) == slice{sliceNum}(1, BotIdx))
                    IdCrit = 1;
                    break;
                end
            end
            if IdCrit == 1
                break;
            end
        end
        if IdCrit == 1
            extraSlice = slice{sliceNum};
            sliceDist(sliceNum - 1) = sliceDist(sliceNum - 1) + zMax - sum(sliceDist);
            slice(sliceNum) = [];
            color(sliceNum) = [];
            sliceNum = sliceNum - 1;
        else
            extraSlice = 0;
            sliceDist(sliceNum - 1) = sliceDist(sliceNum - 1) + zMax - sum(sliceDist);
            slice{sliceNum}(5, : ) = zStd(1);
            slice{1} = [slice{1}, slice{sliceNum}];
            color{1} = [color{1}, color{sliceNum}];
            slice(sliceNum) = [];
            color(sliceNum) = [];
            sliceNum = sliceNum - 1;
        end
    end
    xmin = min(stdLatt(3, : )) - 1;
    xmax = max(stdLatt(3, : )) + 1;
    ymin = min(stdLatt(4, : )) - 1;
    ymax = max(stdLatt(4, : )) + 1;
    if YN == 1
        for SliceIdx = 1 : sliceNum
            ColorList = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'];
            figure;
            hold on;
            for AtomIdx = 1 : size(slice{SliceIdx}, 2)
                scatter(slice{SliceIdx}(3, AtomIdx), slice{SliceIdx}(4, AtomIdx), 'o', ColorList(color{SliceIdx}(AtomIdx)));
            end
%             drawnow;
            axis([xmin xmax ymin ymax]);
            axis equal;
            title(['slice ', num2str(SliceIdx)]);
        end
    end
end

end

