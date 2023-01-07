function [sliceIndices] = DepthsToSliceIndices(sliceDists, stackNum, depths)
%DepthsToSliceIndices.m converts the transmission depths to slice indices.
%   sliceDists -- thickness of slice i or its distance to the next slice or
%       the exit surface;
%   stackNum -- stacking number;
%   depths -- transmission depths;

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

sliceDists = reshape(sliceDists, 1, []);
stackedSliceDists = repmat(sliceDists, 1, stackNum);
stackedSliceDists = cumsum(stackedSliceDists);

sliceIndices = zeros(size(depths));
for depthIdx = 1 : length(depths)
    depthAbsDiff = abs(stackedSliceDists - depths(depthIdx));
    [~, sliceIdx] = min(depthAbsDiff);
    sliceIndices(depthIdx) = sliceIdx;
end

end



