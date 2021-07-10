function [projPotMat, sliceDist] = NanoDiffAppProjPotMat(crystalMatrix,...
    pixelSize, dimPixelNum)
%NanoDiffAppProjPotMat() computes the projected potential matrix as an
%external API for NanoDiff app.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   pixelSize -- sampling pixel size;
%   dimPixelNum -- sampling pixel number per dimension;
% Output:
%   projPot -- projected potential matrix, size: dimPixelNum by dimPixelNum
%       by sliceNum;
%   sliceDist -- slice distance list;

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

if ~isempty(crystalMatrix)
    maxSliceSpacing = 2.0;
    zMax = max(crystalMatrix(5, :)) - min(crystalMatrix(5, :));
    specimenType = 1;
    plotYN = 0;
    [slice, sliceDist, ~] = CrystalSlicing_X(crystalMatrix,...
        crystalMatrix, maxSliceSpacing, zMax, specimenType, plotYN);

    sideLength = pixelSize * dimPixelNum;
    expanNum = [1, 1];
    lattConst = sideLength * [1, 1];

    sliceNum = length(sliceDist);
    projPotMat = zeros(dimPixelNum, dimPixelNum, sliceNum);
    wbHandle = waitbar(0, ['0 / ', num2str(sliceNum)]);
    for sliceIdx = 1 : sliceNum
        tmpSlice = slice{1, sliceIdx};
        tmpSlice(3, :) = tmpSlice(3, :) / sideLength + 0.5;
        tmpSlice(4, :) = tmpSlice(4, :) / sideLength + 0.5;
        projPotMat(:, :, sliceIdx) = MultiProjPot_conv_X(tmpSlice, expanNum,...
            lattConst, sideLength, sideLength, dimPixelNum, dimPixelNum, 1.0e-4);

        waitbar(sliceIdx / sliceNum, wbHandle,...
            [num2str(sliceIdx), ' / ', num2str(sliceNum)]);
    end
    close(wbHandle);
else
    msgbox('Please load xyz file first!');
end

end

