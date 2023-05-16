% prototype_ReduceUnitCell.m
clc;
clear;
close all;
%% init transform
cifFilename = 'tests/advisor/Si_mp-149_symmetrized.cif';
uvw = [1, 1, 1];
[fracCoords, convMat, cutoff] = CrystalAdvisor(cifFilename, uvw);
fullFracCoords = AddEquivAtomSites(fracCoords);

figure;
PlotUnitCell3D(convMat, fullFracCoords);
title('Transformed unit cell (3D)');

%% reduce
maxReduceNum = 4;
% reduce along vec a1
periodicUnit = fullFracCoords;
periodX = 1;
for reduceNum = 2 : maxReduceNum
    glideX = 1 / reduceNum;
    glideUnit = fullFracCoords(:, fullFracCoords(3, :) <= glideX);
    if ~isempty(glideUnit)
        isPeriodic = true;
        for glideIdx = 1 : reduceNum - 1
            tmpGlidedUnit = glideUnit;
            tmpGlidedUnit(3, :) = tmpGlidedUnit(3, :) + glideIdx * glideX;
            lia = ismembertol(tmpGlidedUnit', fullFracCoords', 'ByRows', true);

            if sum(lia, 'all') ~= size(glideUnit, 2)
                isPeriodic = false;
                break;
            end
        end

        if isPeriodic
            periodicUnit = glideUnit;
            periodX = reduceNum;
        end
    end
end
