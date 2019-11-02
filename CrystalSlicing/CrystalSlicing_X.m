function [slice, SliceDist, ExtraSlice] = CrystalSlicing_X(StdLatt, ThermoLatt, DistError, zMax, YN, PlotColor)
%CrystalSlicing_X.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   StdLatt -- standard lattice matrix (without thermo-displacements), uses
%       the most recently accepted crystal info matrix format: [T; P; X; Y;
%       Z], T: atomic types (atomic numbers); P: elemental proportion; X,
%       Y and Z: orthogonal atomic coordinates.
%   ThermoLatt -- lattice with thermo-displacements, adopt the same format
%       of StdLatt;
%   DistError -- the largest error distance to judge whether atoms of
%       different heights be rearranged to one slice;
%   YN -- whether to show each slice: 1 --yes, 0 --no.
%   PlotColor -- colors for each type of atom;
%   NOTE: X denotes an experimental version!
%   Current target: lattice with thermo-displacements.

[zStd, zOrder] = sort(StdLatt(5, : ));
StdLatt = StdLatt( : , zOrder);
ThermoLatt = ThermoLatt( : , zOrder);
PlotColor = PlotColor(zOrder);

SliceNum = 1;
slice{SliceNum} = ThermoLatt( : , 1);
color{SliceNum} = PlotColor(1);
TempZ = zStd(1);
for zIdx = 2 : length(zStd)
    if abs(zStd(zIdx) - TempZ) >= DistError
        SliceDist(SliceNum) = abs(zStd(zIdx) - TempZ);
        SliceNum = SliceNum + 1;
        TempZ = zStd(zIdx);
        ThermoLatt(5, zIdx) = TempZ;
        slice{SliceNum} = ThermoLatt( : , zIdx);
        color{SliceNum} = PlotColor(zIdx);
    else
        ThermoLatt(5, zIdx) = TempZ;
        slice{SliceNum} = [slice{SliceNum}, ThermoLatt( : , zIdx)];
        color{SliceNum} = [color{SliceNum}, PlotColor(zIdx)];
    end
end
% determine whether the last slice should be preserved, deleted or counted
% as the first slice:
if zMax - sum(SliceDist) >= DistError
    SliceDist(SliceNum) = zMax - sum(SliceDist);
    ExtraSlice = 0;
else
    % rearrange the first slice and last slice by ascending order of x,
    % then y, so that whether the last slice is identical to the first
    % slice could be determined
    IdMinDist = 1e-2;
    IdCrit = 0;
    for TopIdx = 1 : size(slice{1}, 2)
        for BotIdx = 1 : size(slice{SliceNum}, 2)
            TempDist = sqrt(sum((slice{1}(3:4, TopIdx) - slice{SliceNum}(3:4, BotIdx)).^2));
            if (TempDist <= IdMinDist) && (slice{1}(1, TopIdx) == slice{4}(1, BotIdx))
                IdCrit = 1;
                break;
            end
        end
        if IdCrit == 1
            break;
        end
    end
    if IdCrit == 1
        ExtraSlice = slice{SliceNum};
        SliceDist(SliceNum - 1) = SliceDist(SliceNum - 1) + zMax - sum(SliceDist);
        slice(SliceNum) = [];
        color(SliceNum) = [];
        SliceNum = SliceNum - 1;
    else
        ExtraSlice = 0;
        SliceDist(SliceNum - 1) = SliceDist(SliceNum - 1) + zMax - sum(SliceDist);
        slice{SliceNum}(5, : ) = zStd(1);
        slice{1} = [slice{1}, slice{SliceNum}];
        color{1} = [color{1}, color{SliceNum}];
        slice(SliceNum) = [];
        color(SliceNum) = [];
        SliceNum = SliceNum - 1;
    end
end
xmin = min(StdLatt(3, : )) - 1;
xmax = max(StdLatt(3, : )) + 1;
ymin = min(StdLatt(4, : )) - 1;
ymax = max(StdLatt(4, : )) + 1;
if YN == 1
    for SliceIdx = 1 : SliceNum
        ColorList = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'];
        figure;
        hold on;
        for AtomIdx = 1 : size(slice{SliceIdx}, 2)
            scatter(slice{SliceIdx}(3, AtomIdx), slice{SliceIdx}(4, AtomIdx), 'o', ColorList(color{SliceIdx}(AtomIdx)));
        end
        axis([xmin xmax ymin ymax]);
        axis equal;
        title(['slice ', num2str(SliceIdx)]);
    end
end

end

