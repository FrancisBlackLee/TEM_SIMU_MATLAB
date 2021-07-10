function [sortCrysMat, sliceInfo, sliceDist] = CrystalSlicing_2(crysMat, distError, YN)
%CrystalSlicing.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   L -- Crystal matrix, where the first row denotes the atomic types, the
%       second row denotes the fractional concentrations and the third to
%       the fifth rows denote the atomic coordinates, whether fractional or
%       orthogonal;
%   DistError -- the largest error distance to judge whether atoms of
%       different heights be rearranged to one slice;
%   YN -- whether to show each slice: 1 --yes, 0 --no.

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

[coordZ, zOrder] = sort(crysMat(5,:));
sortCrysMat = crysMat(:,zOrder);
sliceInfo = 1;
n = 1;
sliceZ = coordZ(1);
for i = 2:length(coordZ)
    if abs(coordZ(i)-sliceZ) >= distError
        sliceInfo = [sliceInfo 1];
        n = n + 1;
        sliceZ = coordZ(i);
    else
        sliceInfo(n) = sliceInfo(n) + 1;
    end
    sortCrysMat(5,i) = sliceZ;
end
sliceDist = zeros(size(sliceInfo));
n = sliceInfo(1);
for i = 1 : length(sliceInfo) - 1
    sliceDist(i) = sortCrysMat(5, n + 1) - sortCrysMat(5, n);
    n = n + sliceInfo(i + 1);
end
sliceDist(i + 1) = sliceDist(1);
% Show the slices
if YN == 1
    n = 1;
    for i = 1:length(sliceInfo)
        figure;
        hold on;
        for j = n:n+sliceInfo(i)-1
            if sortCrysMat(1,j)~=0
                scatter(sortCrysMat(3,j),sortCrysMat(4,j),'o','b');
%                 text(Lp(1,j)+0.2,Lp(2,j),num2str(Lp(4,j)));
            end
        end
        axis([min(sortCrysMat(3,:)) max(sortCrysMat(3,:)) min(sortCrysMat(4,:)) max(sortCrysMat(4,:))]);
        axis equal;
        title(['z= ' num2str(coordZ(n))]);
        n=n+sliceInfo(i);
    end
end

end

