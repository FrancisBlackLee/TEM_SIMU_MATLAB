function [seObjFunc] = SecondaryObjectFunction(keV, fracTypeCoord, expanNum,...
    lattConst, Lx, Ly, Nx, Ny)
%SecondaryObjectFunction.m calculates the secondary electron object
%function for multiple atom. Ref: Wu L, Egerton R F, Zhu Y. Image
%simulation for atomic resolution secondary electron image[J].
%Ultramicroscopy, 2012, 123: 66-73.
%   keV -- energy of incident beam (in keV);
%   fracTypeCoord -- matrix that contains the lattice information, where
%       the first row denotes the atomic types, the second row denotes the
%       elemental proportion, and the third to the fifth row denotes the
%       fractional atomic coordinates. Syntax: [T; P; fracX; fracY; fracZ];
%       Since all the atoms included in this matrix are all on one slice,
%       the fifth row is not required strictly;
%   expanNum -- expansion of the unit cell, syntax: [numX, numY];
%   lattConst -- lattice constants, syntax: [a, b];
%   Lx, Ly, Nx, Ny -- sampling parameters;

% NOTE: the atom removal and lattice expansion are not
% computation-effective enough, which should be optimized later.

% Remove the periodically repeated atoms:
fracTypeCoord = RmvSlcDplAtom_0(fracTypeCoord, 1.0e-8);

% lattice needs to be repeated once more outside the specified expansion
% range, towards both direction (positive & negative) along each dimension,
% to avoid missing atoms.
halfLx = Lx / 2;
halfLy = Ly / 2;
dx = Lx / Nx;
dy = Ly / Ny;
atomNum = size(fracTypeCoord, 2);
seObjFunc = zeros(Ny, Nx);
for atomIdx = 1 : atomNum
    for xExpanIdx = -1 : expanNum(1) + 1
        for yExpanIdx = -1 : expanNum(2) + 1
            atomX = lattConst(1) * (fracTypeCoord(3, atomIdx) + xExpanIdx);
            atomY = lattConst(2) * (fracTypeCoord(4, atomIdx) + yExpanIdx);
            
            % move the center of the expanded lattice to the origin:
            atomX = atomX - halfLx;
            atomY = atomY - halfLy;
            
            tmpSecAmp = SecondaryAmplitude(keV, fracTypeCoord(1, atomIdx),...
                0, 0, Lx, Ly, Nx, Ny);
            tmpSecAmp = imtranslate(tmpSecAmp, [atomX / dx, atomY / dy], 'nearest');
            seObjFunc = seObjFunc + fracTypeCoord(2, atomIdx) * tmpSecAmp.^2;
        end
    end
end

end

