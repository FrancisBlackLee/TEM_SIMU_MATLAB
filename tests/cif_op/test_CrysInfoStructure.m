% test_CrysInfoStructure.m
clc;
clear;
close all;
%% main:
for i = 1 : 5
    strCellArray = {num2str(i), num2str(i + 1)};
    crysInfo.valuedProperty{i, 1} = strCellArray{1};
    crysInfo.valuedProperty{i, 2} = strCellArray{2};
end

for i = 1 : 5
    strCellArray = num2str(i);
    crysInfo.loopProperty{i, 1} = strCellArray;
end

data = {'1', '2', '3', '4', '5'};

for i = 1 : 5
    crysInfo.loopProperty{i, 2} = data{i};
end