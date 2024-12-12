function [dwfs] = MatdynModesToDebyeWallerFactor(qs, bands, eigenVecs, mass, T, thr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 5
    thr = 1e-1;
end

cartAdps = MatdynModesToCartADPs(qs, bands, eigenVecs, mass, T, thr);

dwfs = 8 * pi^2 * cartAdps;

end