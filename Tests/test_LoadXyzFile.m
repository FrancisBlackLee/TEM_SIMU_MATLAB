% test_LoadXyzFile.m
clc;
clear;
close all;
%% main:
filename = 'Tests\Nanocluster.xyz';
filePtr = fopen(filename, 'r');

atomNum = fscanf(filePtr, '%d', 1);
atomTypeList = zeros(1, atomNum);
atomXyzList = zeros(3, atomNum);
for atomIdx = 1 : atomNum
    tmpTypeStr = fscanf(filePtr, '%s', 1);
    atomTypeList(atomIdx) = AtomTypeStrToIdx(tmpTypeStr);
    atomXyzList(1, atomIdx) = fscanf(filePtr, '%f', 1);
    atomXyzList(2, atomIdx) = fscanf(filePtr, '%f', 1);
    atomXyzList(3, atomIdx) = fscanf(filePtr, '%f', 1);
end

fclose(filePtr);

figure;
subplot(2, 2, 1);
scatter3(atomXyzList(1, :), atomXyzList(2, :), atomXyzList(3, :), '.');
title('3D');
axis square;

subplot(2, 2, 2);
scatter(atomXyzList(1, :), atomXyzList(2, :), '.');
title('x-y');
axis square;

subplot(2, 2, 3);
scatter(atomXyzList(1, :), atomXyzList(3, :), '.');
title('x-z');
axis square;

subplot(2, 2, 4);
scatter(atomXyzList(2, :), atomXyzList(3, :), '.');
title('y-z');
axis square;

%% Test:
[testTypeList, testXyzList] = ReadCrystalMakerXyz(filename);
if isequal(testTypeList, atomTypeList) && isequal(testXyzList, atomXyzList)
    msgbox('PASSED');
else
    msgbox('FAILED');
end