% test_CifLineParser.m
clc;
clear;
close all;
%% main:
filename = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_computed.cif'];

fileID = fopen(filename, 'r');

textLine = fgetl(fileID);
while ischar(textLine)
    disp(textLine);
    [strCellArray, lineType, canDelete] = CifLineParser(textLine);
    textLine = fgetl(fileID);
end

fclose(fileID);