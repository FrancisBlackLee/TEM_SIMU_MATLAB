function [adp] = AnisoUToADP(anisoU)
%Convert anisotropic U values to anisotropic displacement parameter tensor.
%   anisoU -- u11, u22, u33, u23, u13, u12

u11 = anisoU(1);
u22 = anisoU(2);
u33 = anisoU(3);
u23 = anisoU(4);
u13 = anisoU(5);
u12 = anisoU(6);
adp = [u11, u12, u13; u12, u22, u23; u13, u23, u33];

end