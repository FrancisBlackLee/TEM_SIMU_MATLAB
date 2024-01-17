function [mappedFxMesh, mappedFyMesh] = MsbtFreqMeshTiltedToPlane(fxMesh, ...
    fyMesh, wavLen, invRotMat)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

fzMesh = sqrt(wavLen^-2 - fxMesh.^2 - fyMesh.^2);
mappedFxMesh = invRotMat(1, 1) * fxMesh + invRotMat(1, 2) * fyMesh +...
    invRotMat(1, 3) * fzMesh;
mappedFyMesh = invRotMat(2, 1) * fxMesh + invRotMat(2, 2) * fyMesh +...
    invRotMat(2, 3) * fzMesh;

end