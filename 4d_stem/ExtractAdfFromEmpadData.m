function [adfImage] = ExtractAdfFromEmpadData(srcFilename, tNum, zNum, adfParams)
%ExtractAdfFromEmpadData.m extracts the annular dark field image from the
%EMPAD data, however, the integration range does not have to be dark field,
%bright field is also allowed.
%   srcFilename -- filename of EMPAD data;
%   tNum, zNum -- number of scanning positions per dimension, t: dimension
%       1 (vertical), z: dimension 2 (horizontal);
%   adfParams:
%       adfParams.innerRadius -- inner radius of the ADF detector (in
%           pixel);
%       adfParams.outerRadius -- outer radius of the ADF detector (in
%           pixel);
%       adfParams.scaling -- scaling of the ADF detector, in order to
%           reduce the loss of accuracy due to large pixel size;
%       adfParams.interpolationMethod -- interpolation method used in
%           scaling, translation, and back-scaling of the ADF detector;

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

[fileID, errMsg] = fopen(srcFilename);
if fileID == -1
    error(errMsg);
end

imageSize = [128, 130];
cropRangeY = 3 : 126;
cropRangeX = 3 : 126;

cropNy = 124;
cropNx = 124;
cropAxisX = -cropNx / 2 : cropNx / 2 - 1;
cropAxisY = -cropNy / 2 : cropNy / 2 - 1;
[cropMeshX, cropMeshY] = meshgrid(cropAxisX, cropAxisY);

adfImage = zeros(tNum, zNum);
wbHandle = waitbar(0, 'processing...');
for tIdx = 1 : tNum
    for zIdx = 1 : zNum
        tmpRawImage = fread(fileID, imageSize, 'float', 'ieee-le');
        tmpCropImage = tmpRawImage(cropRangeY, cropRangeX);
        [tmpCenterX, tmpCenterY] = FindDiskCenter(tmpCropImage);
        adfDetector = GenerateAdfDetector(tmpCenterX, tmpCenterY);
        adfImage(tIdx, zIdx) = sum(tmpCropImage .* adfDetector, 'all');
    end
    wbMsg = [num2str(tIdx), ' / ', num2str(tNum), ' completed'];
    waitbar(tIdx / tNum, wbHandle, wbMsg);
end
delete(wbHandle);

fclose(fileID);


% nested functions
    function [centerX, centerY] = FindDiskCenter(I)
        normI = I / sum(I, 'all');
        centerX = sum(normI .* cropMeshX, 'all');
        centerY = sum(normI .* cropMeshY, 'all');
    end

    function adfDetector = GenerateAdfDetector(centerX, centerY)
        scaling = adfParams.scaling;
        
        scaledNx = scaling * cropNx;
        scaledNy = scaling * cropNy;
        scaledAxisX = -scaledNx / 2 : scaledNx / 2 - 1;
        scaledAxisY = -scaledNy / 2 : scaledNy / 2 - 1;
        [scaledMeshX, scaledMeshY] = meshgrid(scaledAxisX, scaledAxisY);
        
        scaledCenterX = scaling * centerX;
        scaledCenterY = scaling * centerY;
        scaledInnerR = scaling * adfParams.innerRadius;
        scaledOuterR = scaling * adfParams.outerRadius;
        scaledMeshR = sqrt((scaledMeshX - scaledCenterX).^2 +...
            (scaledMeshY - scaledCenterY).^2);
        scaledAdfDetector = 1.0 * ((scaledMeshR > scaledInnerR) &...
            (scaledMeshR < scaledOuterR));
        adfDetector = imresize(scaledAdfDetector, 1 / scaling,...
            adfParams.interpolationMethod);
    end

end

