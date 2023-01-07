function [fracInten] = VtemlabPacbedQstem(pacbeds, detector)
%VtemlabPacbedQstem.m analyses the pacbed matrix for quantitative scanning
%transmission electron microscopy.
% Input:
%   pacbeds -- position averaged convergent beam electron diffraction
%       patterns, formed as a 3D matrix, whose third dimension denotes the
%       number of patterns;
%   detector -- integral detector, i.e., HAADF detector;
% Output:
%   fracInten -- fractional intensities.

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

imageNum = size(pacbeds, 3);
fracInten = zeros(1, imageNum);
for imageIdx = 1 : imageNum
    tmpImage = pacbeds(:, :, imageIdx);
    fracInten(imageIdx) = sum(tmpImage .* detector, 'all') / sum(tmpImage, 'all');
end

end


