function [a1p, a2p, rotMat, invRotMat] = NonorthoBasesReorient(a1, a2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a1 = reshape(a1, [], 1);
a2 = reshape(a2, [], 1);
a3 = cross(a1, a2);
cosA1A2 = dot(a1, a2) / (norm(a1) * norm(a2));
sinA1A2 = sqrt(1 - cosA1A2^2);

a1p = [norm(a1); 0; 0];
a2p = [norm(a2) * cosA1A2; norm(a2) * sinA1A2; 0];
a3p = cross(a1p, a2p);

rotMat = [a1p, a2p, a3p] / [a1, a2, a3];
invRotMat = inv(rotMat);

end