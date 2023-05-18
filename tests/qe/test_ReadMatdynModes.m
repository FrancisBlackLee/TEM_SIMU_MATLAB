% test_ReadMatdynModes.m
clc;
clear;
close all;
%% main:
filename = 'tests/qe/matdyn.modes';
[qs, bands, eigenVecs] = ReadMatdynModes(filename);
