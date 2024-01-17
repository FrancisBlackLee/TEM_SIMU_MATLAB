function [detector] = MeshedAnnularDetector(fxMesh, fyMesh, wavLen, ...
    innerAngle, outerAngle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

freqSqr = fxMesh.^2 + fyMesh.^2;
detector = (freqSqr >= (innerAngle * 1.0e-3 / wavLen)^2) &...
    (freqSqr <= (outerAngle * 1.0e-3 / wavLen)^2);

end