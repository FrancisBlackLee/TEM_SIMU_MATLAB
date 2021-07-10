function [slice, sliceDist] = CrystalSlicing_3(crysMat, distError, zMax, YN, plotColor)
%CrystalSlicing.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   crysMat -- Crystal matrix, where the first row denotes the atomic types, the
%       second row denotes the fractional concentrations and the third to
%       the fifth rows denote the atomic coordinates, whether fractional or
%       orthogonal;
%   distError -- the largest error distance to judge whether atoms of
%       different heights be rearranged to one slice;
%   YN -- whether to show each slice: 1 --yes, 0 --no.
%   NOTE: this version was derived from the previous experimental version,
%   designed for LPCMO slicing.

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

atomNum = size(crysMat, 2);
[coordZ, zOrder] = sort(crysMat(5,:));
sortCrysMat = crysMat(:,zOrder);
plotColor = plotColor(zOrder);
sliceInfo = 1;
n = 1;
slice{n} = sortCrysMat( : , 1);
sliceZ = coordZ(1);
for i = 2:length(coordZ)
    if abs(coordZ(i)-sliceZ) >= distError
        sliceInfo = [sliceInfo 1];
        n = n + 1;
        slice{n} = sortCrysMat( : , i);
        sliceZ = coordZ(i);
    else
        sliceInfo(n) = sliceInfo(n) + 1;
        sortCrysMat(5, i) = sliceZ;
        slice{n} = [slice{n}, sortCrysMat( : , i)];
    end
end
sliceDist = zeros(size(sliceInfo));
n = sliceInfo(1);
for i = 1 : length(sliceInfo) - 1
    sliceDist(i) = sortCrysMat(5, n + 1) - sortCrysMat(5, n);
    n = n + sliceInfo(i + 1);
end
if zMax - sum(sliceDist(1 : i)) >= distError
    sliceDist(i + 1) = zMax - sum(sliceDist(1 : i));
else
    sliceDist(i) = sliceDist(i) + zMax - sum(sliceDist(1 : i));
    slice{i+1}(5, : ) = coordZ(1);
    sortCrysMat(5, atomNum - sliceInfo(i + 1) + 1 : atomNum) = sortCrysMat(5, 1);
    sortCrysMat = [sortCrysMat( : , atomNum - sliceInfo(i + 1) + 1 : atomNum), sortCrysMat( : , 1 : atomNum - sliceInfo(i + 1))];
    plotColor = [plotColor(atomNum - sliceInfo(i + 1) + 1 : atomNum), plotColor(1 : atomNum - sliceInfo(i + 1))];
    slice{1} = [slice{1}, slice{i+1}];
    slice(i + 1) = [];
    sliceDist(i + 1) = [];
    sliceInfo(1) = sliceInfo(1) + sliceInfo(i + 1);
    sliceInfo(i + 1) = [];
end
% Show the slices
if YN == 1
    n = 1;
    for i = 1:length(sliceInfo)
        figure;
        hold on;
        for j = n:n+sliceInfo(i)-1
            if sortCrysMat(1,j)~=0
                % No more than 8 types of color
                Colors = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'];
                scatter(sortCrysMat(3, j), sortCrysMat(4, j), 'o', Colors(plotColor(j)));
            end
        end
        axis([min(sortCrysMat(3,:)) max(sortCrysMat(3,:)) min(sortCrysMat(4,:)) max(sortCrysMat(4,:))]);
        axis equal;
        title(['z= ' num2str(coordZ(n))]);
        n=n+sliceInfo(i);
    end
end

end

