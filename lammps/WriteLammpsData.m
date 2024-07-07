function WriteLammpsData(filename, convMat, axyz, thr)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

if nargin == 3
    thr = 1.0e-6;
end

% if it is triclinic, use the general form
if convMat(2, 1) > thr || convMat(3, 1) > thr || convMat(3, 2) > thr
    error("The input conversion matrix is not the general form of a triclinic system");
end

fileId = fopen(filename, "w");
if fileId == -1
    error("Failed to open file");
else
    lx = convMat(1, 1);
    xy = convMat(1, 2);
    xz = convMat(1, 3);
    ly = convMat(2, 2);
    yz = convMat(2, 3);
    lz = convMat(3, 3);

    atypes = axyz(1, :);
    nAtom = length(atypes);
end

end