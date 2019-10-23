function [NewMat] = DelRepAtom_0(CrysInfoMat)
%DelRepAtom_0.m deletes the repeated atoms in a crystal information matrix
%exported by CrystalMaker.
%   CrysInfoMat -- an loaded matrix containing the crystal positional
%       information, exported by CrystalMaker and later processed to
%       include only numerical data, whose syntax is:
%       [T1,     ..., TN;
%        T1,     ..., TN;
%        xfrac1, ..., xfracN;
%        yfrac1, ..., yfracN;
%        zfrac1, ..., zfracN;
%        xor1,   ..., xorN;
%        yor1,   ..., yorN;
%        zor1,   ..., zorN];
%       T denotes atomic type represented by its atomic number;
%       xfrac, yfrac, and zfrac denote the fractional coordinates;
%       xor, yor, and zor denote the orthogonal coordinates;
%   Note: because the precision of the exported floating-point number from
%   CrystalMaker is 1e-6, so the error distance must be no greater than
%   1e-6

NewMat = CrysInfoMat;
% search xfrac and yfrac periodic duplication:
RefIdx = 1; % Current reference index
while RefIdx < size(NewMat, 2)
    CurIdx = RefIdx + 1; % Current checking index
    while CurIdx <= size(NewMat, 2)
        if (NewMat(1, CurIdx) == NewMat(1, RefIdx)) && (((abs(abs(NewMat(3, CurIdx) - NewMat(3, RefIdx)) - 1) <= 1e-7)...
                && (abs(NewMat(4, CurIdx) - NewMat(4, RefIdx)) <= 1e-7)) || ((abs(NewMat(3, CurIdx) - NewMat(3, RefIdx)) <= 1e-7)...
                && (abs(abs(NewMat(4, CurIdx) - NewMat(4, RefIdx)) - 1) <= 1e-7))) && (abs(NewMat(5, CurIdx) - NewMat(5, RefIdx)) <= 1e-7)
            NewMat( : , CurIdx) = [];
        else
            CurIdx = CurIdx + 1;
        end
    end
    RefIdx = RefIdx + 1;
end

end

