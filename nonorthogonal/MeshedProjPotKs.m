function [projPot] = MeshedProjPotKs(fxMesh, fyMesh, stretchCoeff, types, ...
    occus, axs, ays)
%MeshedProjPotKs calculates projected potential of atoms in reciprocal space
%using a pair of predefined fx and fy meshes.
% [projPot] = MeshedProjPotKs(fxMesh, fyMesh, stretchCoeff, types, occus, axs, ays)
%   fxMesh -- mesh of fx coordinates;
%   fyMesh -- mesh of fy coordinates;
%   stretchCoeff -- 
%   types -- atomic types;
%   occus -- atomic occupations;
%   axs -- atomic x coordinates;
%   ays -- atomic y coordinates;
%   projPot -- projected potential;

a0 = 0.529;
e0 = 14.4;
scaleCoeff = 2.0 * pi * e0 * a0;

sf = zeros(size(fxMesh));
frMesh = sqrt(fxMesh.^2 + fyMesh.^2);

na = length(types);
for ia = 1 : na
    tmpSf = scaleCoeff * ScatteringFactor(types(ia), frMesh);
    tmpKer = occus(ia) * exp(-1i * 2.0 * pi * (fxMesh * axs(ia) + ...
        fyMesh * ays(ia)));
    sf = sf + tmpKer .* tmpSf;
end

projPot = stretchCoeff * real(ifftshift(ifft2(fftshift(sf))));

end