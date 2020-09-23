function [newAtomSiteMat] = AddEquivAtomSites(atomSiteMat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

initAtomNum = size(atomSiteMat, 2);

newAtomSiteMat = repmat(atomSiteMat, 1, 27);
for xGlide = -1 : 1
    start_x = (xGlide + 1) * initAtomNum * 9;
    for yGlide = -1 : 1
        start_y = (yGlide + 1) * initAtomNum * 3;
        for zGlide = -1 : 1
            start_z = (zGlide + 1) * initAtomNum;
            startIdx = start_x + start_y + start_z + 1;
            range = startIdx : startIdx + initAtomNum - 1;
            newAtomSiteMat(3, range) = newAtomSiteMat(3, range) + xGlide;
            newAtomSiteMat(4, range) = newAtomSiteMat(4, range) + yGlide;
            newAtomSiteMat(5, range) = newAtomSiteMat(5, range) + zGlide;
        end
    end
end

tolerance = 1e-8;

index = find(newAtomSiteMat(3, :) > 1 + tolerance);
newAtomSiteMat(3, index) = mod(newAtomSiteMat(3, index), 1);
index = find(newAtomSiteMat(3, :) < 0 - tolerance);
newAtomSiteMat(3, index) = newAtomSiteMat(3, index) - floor(newAtomSiteMat(3, index));

index = find(newAtomSiteMat(4, :) > 1 + tolerance);
newAtomSiteMat(4, index) = mod(newAtomSiteMat(4, index), 1);
index = find(newAtomSiteMat(4, :) < 0 - tolerance);
newAtomSiteMat(4, index) = newAtomSiteMat(4, index) - floor(newAtomSiteMat(4, index));

index = find(newAtomSiteMat(5, :) > 1 + tolerance);
newAtomSiteMat(5, index) = mod(newAtomSiteMat(5, index), 1);
index = find(newAtomSiteMat(5, :) < 0 - tolerance);
newAtomSiteMat(5, index) = newAtomSiteMat(5, index) - floor(newAtomSiteMat(5, index));

newAtomSiteMat = newAtomSiteMat';
newAtomSiteMat = uniquetol(newAtomSiteMat, tolerance, 'ByRows', true');
newAtomSiteMat = newAtomSiteMat';

end

