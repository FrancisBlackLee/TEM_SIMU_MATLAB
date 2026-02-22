function [seObjFunc] = SecondaryObjectFunctionCart_gpu(keV, cartTypeCoord, ...
    Lx, Ly, Nx, Ny)

%SecondaryObjectFunctionCart.m calculates the secondary electron object
%function for multiple atom. Ref: Wu L, Egerton R F, Zhu Y. Image
%simulation for atomic resolution secondary electron image[J].
%Ultramicroscopy, 2012, 123: 66-73.
%   keV -- energy of incident beam (in keV);
%   cartTypeCoord -- matrix that contains the lattice information, where
%       the first row denotes the atomic types, the second row denotes the
%       elemental proportion, and the third to the fifth row denotes the
%       cartesian atomic coordinates. Syntax: [T; P; X; Y; Z];
%       Since all the atoms included in this matrix are all on one slice,
%       the fifth row is not required strictly;
%   Lx, Ly, Nx, Ny -- sampling parameters;

% NOTE: the atom removal and lattice expansion are not
% computation-effective enough, which should be optimized later.

nAtom = size(cartTypeCoord, 2);
seObjFunc = zeros(Ny, Nx, "single", "gpuArray");
for iAtom = 1 : nAtom
    secAmp = SecondaryAmplitude_gpu(keV, cartTypeCoord(1, iAtom), ...
        cartTypeCoord(3, iAtom), cartTypeCoord(4, iAtom), Lx, Ly, Nx, Ny);
    seObjFunc = seObjFunc + cartTypeCoord(2, iAtom) * secAmp.^2;
end

end