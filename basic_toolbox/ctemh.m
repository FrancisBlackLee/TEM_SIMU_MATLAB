function y = ctemh(k, params)
% I just copied it from E. J. Kirkland's book, I was lazy...
%
% MATLAB function ctemh.m to calculate CTEM bright
% field phase contrast transfer function with partial
% coherence for weak phase objects
%
% input array k has the spatial freq. values (in 1/A)
%
%
%
%
%
%
% input array params has the optical parameters
% params = [Cs, df, kev, ddf, beta];
% input type = 0 for phase contrast
% and 1 for amplitude contrast
% output array contains the transfer function vs k
% params.Cs3,5 = spherical aberration (in mm)
% params.df = defocus (in Angstroms)
% params.kev = electron energy (in keV)
% params.ddf = chrom. aberr. def. spread (in Angst.)
% params.beta = spread in illum. angles (in mrad)
%
% reference
% R. H. Wade and J. Frank, Optik 49 (1977) p.81
%
Cs3 = params.Cs3*1.0e7;
Cs5 = params.Cs5*1.0e7;
df = params.df;
kev = params.kev;
ddf = params.ddf;
beta = params.beta*0.001;
mo = 511.0;
hc = 12.3986;
wav = (2*mo)+kev;
% electron rest mass in keV
% in keV-Angstroms
wav = hc/sqrt(wav*kev);
wavsq = wav*wav;
w1 = pi*Cs3*wavsq*wav;
w2 = pi*wav*df;
w3 = pi*Cs5*wavsq*wavsq*wav;
e0 = (pi*beta*ddf)^2;
k2 = k .* k;
wr = ((w3.*k2+w1).*k2-w2).*k*beta/wav;
wi = pi*wav*ddf.*k2;
wi = wr.*wr + 0.25.*wi.*wi;
wi = exp(-wi./(1+e0.*k2));
wr = w3*(1-2.0*e0.*k2)/3.0;
wr = wr.*k2 + 0.5*w1.*(1-e0.*k2);
wr = (wr.*k2- w2).*k2./(1+e0.*k2);
y = exp(1i * wr) .* wi;

end