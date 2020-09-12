% test_SelectFile.m
clc;
clear;
close all;
%% main:
[file, path] = uigetfile('*.xyz');
if isequal(file, 0)
    msgbox('User selected cancel');
else
    msgbox(['User selected ', fullfile(path, file)]);
end