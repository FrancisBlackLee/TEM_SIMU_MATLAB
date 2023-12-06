% prototype_extract_unitcell_image_1.m
clc;
clear;
close all;
%% load raw data:
rawRoi = load("tests/geometry/raw1_roi.txt");
rawPos = load("tests/geometry/roi_pos.txt");

[rawNy, rawNx] = size(rawRoi);
rawX = InitAxis(rawNx, rawNx);
rawY = InitAxis(rawNy, rawNy);
[rawXMesh, rawYMesh] = meshgrid(rawX, rawY);

rawPos = rawPos + 1;
grayRawRoi = mat2gray(rawRoi);
markedRawRoi = insertMarker(grayRawRoi, rawPos, "circle", "Size", 8, "Color", 'red');
figure;
imshow(markedRawRoi);
title("marked raw roi");

h = drawpolygon("FaceAlpha", 0);
unitcellPos = h.Position;
h.Color = 'yellow';

vec1 = unitcellPos(2, :) - unitcellPos(1, :);
vec2 = unitcellPos(3, :) - unitcellPos(1, :);
vec3 = unitcellPos(4, :) - unitcellPos(1, :);

%% load unit cell:
[unitcell, convMat, cutoff] = CrystalAdvisor("tests/geometry/MoS2_mp-2815_symmetrized.cif", ...
    [0, 0, 1], true);

unitcell = GlideUnitCell(unitcell, [1/3, 0, 0]);
convMat = rotz(90) * convMat;

figure;
PlotUnitCell2D(convMat, unitcell);

a1 = convMat(1 : 2, 1)';
a2 = convMat(1 : 2, 2)';
n1 = round(norm(a1) / 0.1) + 1;
n2 = round(norm(a2) / 0.1) + 1;
gridA1 = linspace(0, 1, n1);
gridA2 = linspace(0, 1, n2);
[meshA1, meshA2] = meshgrid(gridA1, gridA2);
xMesh = a1(1) * meshA1 + a2(1) * meshA2;
yMesh = a1(2) * meshA1 + a2(2) * meshA2;

refPos = [0, 0; a1; a2; a1 + a2] / 0.1;

%% find all unit cells
[permRefIdx, sortedRefPos] = matchpoints(unitcellPos, refPos);
nPos = size(rawPos, 1);
thr = 5;
avRoi = zeros(size(rawRoi));
avNum = 0;
for iPos = 1 : nPos
    tmpPos0 = rawPos(iPos, :);
    % find tmpPos0 + vec1
    tmpPos1 = tmpPos0 + vec1;
    tmpDiff = rawPos - tmpPos1;
    tmpDists = sqrt(tmpDiff(:, 1).^2 + tmpDiff(:, 2).^2);
    [mr, mi] = min(tmpDists);
    if mr > thr
        continue;
    end

    tmpPos1 = rawPos(mi, :);

    % find tmpPos0 + vec2
    tmpPos2 = tmpPos0 + vec2;
    tmpDiff = rawPos - tmpPos2;
    tmpDists = sqrt(tmpDiff(:, 1).^2 + tmpDiff(:, 2).^2);
    [mr, mi] = min(tmpDists);
    if mr > thr
        continue;
    end

    tmpPos2 = rawPos(mi, :);

    % find tmpPos0 + vec3
    tmpPos3 = tmpPos0 + vec3;
    tmpDiff = rawPos - tmpPos3;
    tmpDists = sqrt(tmpDiff(:, 1).^2 + tmpDiff(:, 2).^2);
    [mr, mi] = min(tmpDists);
    if mr > thr
        continue;
    end

    tmpPos3 = rawPos(mi, :);

    % moving points
    mvPos = [tmpPos0; tmpPos1; tmpPos2; tmpPos3];
    tform = fitgeotrans(mvPos, sortedRefPos, 'affine');
    tmpRegRoi = imwarp(rawRoi, tform, 'OutputView', imref2d(size(markedRawRoi)), 'interp', 'linear');
    avRoi = avRoi + tmpRegRoi;
    avNum = avNum + 1;
end
avRoi = avRoi / avNum;


