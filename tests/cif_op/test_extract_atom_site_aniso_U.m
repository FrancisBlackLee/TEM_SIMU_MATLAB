% test_extract_atom_site_aniso_U.m
clc;
clear;
close all;
%% main:
crysInfo = LoadCif('tests/cif_op/tdispmat.cif');
[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
convMat = ConversionMatrix(cellLengths, cellAngles);
atomSites = ExtractAtomSiteFromCrysInfo(crysInfo);
anisoUs = ExtractAtomSiteAnisoUFromCrysInfo(crysInfo);