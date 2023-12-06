function [unitCell, convMat] = CrysInfoToUnitCell(crysInfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
unitCell = ExtractAtomSiteFromCrysInfo(crysInfo);
convMat = ConversionMatrix(cellLengths, cellAngles);

end