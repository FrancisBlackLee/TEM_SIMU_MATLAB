function [avRoi] = AverageUnitCellImage_X(rawRoi, rawPos, scanNx, scanNy, pixThr, stdThr, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 4
    pixThr = 10;
    stdThr = 0.15;
end

% draw unit cell shape
grayRawRoi = mat2gray(rawRoi);
markedRawRoi = insertMarker(grayRawRoi, rawPos, varargin{:});
figure;
imshow(markedRawRoi);
title("marked raw roi");

h = drawpolygon("FaceAlpha", 0);
unitcellPos = h.Position;
h.Color = 'yellow';

vec1 = unitcellPos(2, :) - unitcellPos(1, :);
vec2 = unitcellPos(3, :) - unitcellPos(1, :);
vec3 = unitcellPos(4, :) - unitcellPos(1, :);

% find all unit cells
refPos = [1, scanNy; scanNx, scanNy; scanNx, 1; 1, 1];
[permRefIdx, sortedRefPos] = matchpoints(unitcellPos, refPos);
nPos = size(rawPos, 1);
avRoi = zeros(scanNy, scanNx);
avNum = 0;
figure;
for iPos = 1 : nPos
    tmpPos0 = rawPos(iPos, :);
    % find tmpPos0 + vec1
    tmpPos1 = tmpPos0 + vec1;
    tmpDiff = rawPos - tmpPos1;
    tmpDists = sqrt(tmpDiff(:, 1).^2 + tmpDiff(:, 2).^2);
    [mr, mi] = min(tmpDists);
    if mr > pixThr
        continue;
    end

    tmpPos1 = rawPos(mi, :);

    % find tmpPos0 + vec2
    tmpPos2 = tmpPos0 + vec2;
    tmpDiff = rawPos - tmpPos2;
    tmpDists = sqrt(tmpDiff(:, 1).^2 + tmpDiff(:, 2).^2);
    [mr, mi] = min(tmpDists);
    if mr > pixThr
        continue;
    end

    tmpPos2 = rawPos(mi, :);

    % find tmpPos0 + vec3
    tmpPos3 = tmpPos0 + vec3;
    tmpDiff = rawPos - tmpPos3;
    tmpDists = sqrt(tmpDiff(:, 1).^2 + tmpDiff(:, 2).^2);
    [mr, mi] = min(tmpDists);
    if mr > pixThr
        continue;
    end

    tmpPos3 = rawPos(mi, :);

    % moving points
    mvPos = [tmpPos0; tmpPos1; tmpPos2; tmpPos3];
    tform = fitgeotrans(mvPos, sortedRefPos, 'affine');
    tmpRegRoi = imwarp(rawRoi, tform, 'OutputView', imref2d(size(markedRawRoi)), 'interp', 'linear');

    % compare the tmpRoi with avRoi
    tmpRoi = tmpRegRoi(1 : scanNy, 1 : scanNx);
    if avNum > 0
        normAvRoi = avRoi / mean(avRoi, 'all');
        normTmpRoi = tmpRoi / mean(tmpRoi, 'all');
        roiDiff = normTmpRoi - normAvRoi;
        roiStd = sqrt(mean(roiDiff.^2, 'all'));
        imshowpair(avRoi, tmpRoi, 'montage');
        title(['std=', num2str(roiStd)]);
        drawnow;
        if roiStd < stdThr
            avRoi = avRoi + tmpRoi;
            avNum = avNum + 1;
        end
    else
        avRoi = avRoi + tmpRoi;
        avNum = avNum + 1;
    end
end
avRoi = avRoi / avNum;

end