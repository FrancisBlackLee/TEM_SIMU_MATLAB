function [mappedFxMesh, mappedFyMesh] = MsbtFreqMeshTiltedToPlane(fxMesh, ...
    fyMesh, wavLen, invRotMat)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

frMesh = sqrt(fxMesh.^2 + fyMesh.^2);
fzMesh = 1 / wavLen * sqrt(1 ./ (1 + wavLen^2 * frMesh.^2));
mappedFxMesh = invRotMat(1, 1) * fxMesh + invRotMat(1, 2) * fyMesh +...
    invRotMat(1, 3) * fzMesh;
mappedFyMesh = invRotMat(2, 1) * fxMesh + invRotMat(2, 2) * fyMesh +...
    invRotMat(2, 3) * fzMesh;

end