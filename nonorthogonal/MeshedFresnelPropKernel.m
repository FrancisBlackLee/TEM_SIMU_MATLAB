function propKernel = MeshedFresnelPropKernel(fxMesh, fyMesh, wavLen, propDist, ...
    tiltX, tiltY)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 4
    propKernel = exp(-1i * pi * wavLen * propDist * (fxMesh.^2 + fyMesh.^2));
elseif nargin == 6
    tx = tiltX * 1.0e-3;
    ty = tiltY * 1.0e-3;
    propKernel = exp(-1i * pi * wavLen * propDist * (fxMesh.^2 + fyMesh.^2) +...
        1i * 2 * pi * propDist * (fxMesh * tan(tx) + fyMesh * tan(ty)));
else
    error('Incorrect number of arguments.');
end

end