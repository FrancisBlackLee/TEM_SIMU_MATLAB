function aperture = MeshedAnnularAperture(fxMesh, fyMesh, wavLen, r0, r1, pol, azi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 5
    freqSqr = fxMesh.^2 + fyMesh.^2;
    aperture = (freqSqr > (r0 * 1.0e-3 / wavLen)^2  &...
        freqSqr < (r1 * 1.0e-3 / wavLen)^2);
elseif nargin == 7
    tiltFreq = pol * 1.0e-3 / wavLen;
    tiltFreqX = tiltFreq * cosd(azi);
    tiltFreqY = tiltFreq * sind(azi);
    freqSqr = (fxMesh - tiltFreqX).^2 + (fyMesh - tiltFreqY).^2;
    aperture = (freqSqr > (r0 * 1.0e-3 / wavLen)^2 &...
        freqSqr < (r1 * 1.0e-3 / wavLen)^2);
else
    error("Incorrect number of input arguments");
end

end