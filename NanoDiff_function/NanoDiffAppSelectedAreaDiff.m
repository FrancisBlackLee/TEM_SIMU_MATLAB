function [diffPatt] = NanoDiffAppSelectedAreaDiff(keV, projPotMat,...
    sliceDist, pixelSize, dimPixelNum)
%NanoDiffAppSelectedAreaDiff() computes the selected area diffraction
%pattern as an external API for NanoDiff app.
% Input:
%   keV -- electron beam energy;
%   projPotMat -- projected potential matrix;
%   sliceDist -- slice distance list;
% Output:
%   diffPatt -- diffraction pattern;

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

sideLength = pixelSize * dimPixelNum;
interCoeff = InteractionCoefficient(keV);
sliceNum = length(sliceDist);

transFuncMat = 1i * ones(dimPixelNum, dimPixelNum, sliceNum);
for sliceIdx = 1 : sliceNum
    transFuncMat(:, :, sliceIdx) =...
        exp(1i * interCoeff * projPotMat(:, :, sliceIdx) / 1e3);
end

wave = ones(dimPixelNum, dimPixelNum);
stackNum = 1;
wave = multislice_X(wave, keV, sideLength, sideLength, transFuncMat,...
    sliceDist, stackNum);

wave = ifftshift(fft2(fftshift(wave))) * pixelSize^2;
diffPatt = abs(wave.^2);
diffPatt = log(1 + 0.1 * diffPatt);

end

