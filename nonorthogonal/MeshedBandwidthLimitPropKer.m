function limited = MeshedBandwidthLimitPropKer(unlimited, fxMesh, fyMesh, bwl)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

frMesh = sqrt(fxMesh.^2 + fyMesh.^2);
limited = (frMesh <= bwl) .* unlimited;

end