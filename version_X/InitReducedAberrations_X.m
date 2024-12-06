function [aberrations] = InitReducedAberrations_X()
%InitReducedAberrations_X inits aberrations including only Cs3, Cs5 and
%defocus.
%   aberrations:
%       aberrations.Cs3 -- 3rd order spherical aberration in mm;
%       aberrations.Cs5 -- 5th order spherical aberration in mm;
%       aberrations.df -- defocus in angstrom (a negative defocus to 
%           eliminate the effect of spherical aberrations, different from 
%           old versions);

aberrations.Cs3 = 0;
aberrations.Cs5 = 0;
aberrations.df = 0;

end