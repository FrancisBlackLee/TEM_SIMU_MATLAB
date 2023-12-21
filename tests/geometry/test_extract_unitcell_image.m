% test_extract_unitcell_image.m
clc;
clear;
close all;
%% main:
rawRoi = load("tests/geometry/raw1_roi.txt");
rawPos = load("tests/geometry/roi_pos.txt");

rawPos = rawPos + 1;

scanNx = 36;
scanNy = 62;
thr = 10;

avRoi = AverageUnitCellImage(rawRoi, rawPos, scanNx, scanNy, thr);

figure;
imshow(avRoi, []);