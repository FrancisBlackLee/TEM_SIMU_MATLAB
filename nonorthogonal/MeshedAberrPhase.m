function phase = MeshedAberrPhase(aberrs, wavLen, fxMesh, fyMesh)
%MeshedAberrPhase computes the aberration phase shift with the preset
%meshes.
% phase = MeshedAberrPhase(aberrs, aberrAngles, wavLen, fxMesh, fyMesh)
%   aberrs -- aberrations:
%       C1  Defocus
%       A1  2-Fold astigmatism
%       B2  Axial coma
%       A2  3-Fold astigmatism
%       C3  3rd order spherical aberration
%       S3  Axial star aberration
%       A3  4-Fold astigmatism
%       B4  4th order axial coma
%       D4  3-Lobe aberration
%       A4  5-Fold astigmatism
%       C5  5th order spherical aberration
%       S5  5th order axial star
%       R5  5th order rosette
%       A5  6-Fold astigmatism
%
%       A1_angle  2-Fold astigmatism
%       B2_angle  Axial coma
%       A2_angle  3-Fold astigmatism
%       S3_angle  Axial star aberration
%       A3_angle  4-Fold astigmatism
%       B4_angle  4th order axial coma
%       D4_angle  3-Lobe aberration
%       A4_angle  5-Fold astigmatism
%       S5_angle  5th order axial star
%       R5_angle  5th order rosette
%       A5_angle  6-Fold astigmatism
%   wavLen -- electron wavelength
%   fxMesh, fyMesh -- preset meshes of spatial frequencies.

a{1} = [aberrs.C1, aberrs.A1];
a{2} = [aberrs.B2, aberrs.A2];
a{3} = [aberrs.C3, aberrs.S3, aberrs.A3];
a{4} = [aberrs.B4, aberrs.D4, aberrs.A4];
a{5} = [aberrs.C5, aberrs.S5, aberrs.R5, aberrs.A5];

r{1} = [0, aberrs.A1_angle];
r{2} = [aberrs.B2_angle, aberrs.A2_angle];
r{3} = [0, aberrs.S3_angle, aberrs.A3_angle];
r{4} = [aberrs.B4_angle, aberrs.D4_angle, aberrs.A4_angle];
r{5} = [0, aberrs.S5_angle, aberrs.R5_angle, aberrs.A5_angle];

aberrNum = [2, 2, 3, 3, 4]; % number of aberrations in each aberration array

polAng = wavLen * sqrt(fxMesh.^2 + fyMesh.^2);
aziAng = atan2(fyMesh, fxMesh);
% Aberr{n}(Idx) = Cnm, m = (Idx - 1) * 2 + (n - 1) % 2
phase = zeros(size(fxMesh));
for n = 1 : 5
    for Idx = 1 : aberrNum(n)
        m = (Idx - 1) * 2 + mod(n - 1, 2);
        phase = phase + a{n}(Idx) * polAng.^(n + 1)...
            .* cos(m * (aziAng - deg2rad(r{n}(Idx)))) / (n + 1);
    end
end
phase = 2 * pi / wavLen * phase;

end