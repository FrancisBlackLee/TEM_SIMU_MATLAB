function [rotMat] = RotationOperator(direction)
%RotationOperator() computes the martix form of the rotation operator given
%the view or projection direction.
% Input:
%   direction -- view or projection direction, vector with 3 elements;
% Output:
%   rotMat -- rotation operator matrix for cartesian coordinates:
%       usage: [xp; yp; zp] = rotMat * [x; y; z];

newVecZ = reshape(direction, 1, []);
newVecX = [0, newVecZ(3), -newVecZ(2)];
if ~any(newVecX)
    newVecX = [0, 1, 0];
    newVecZ = newVecZ / norm(newVecZ);
else
    % normalization:
    newVecZ = newVecZ / norm(newVecZ);
    newVecX = newVecX / norm(newVecX);
end

newVecY = cross(newVecZ, newVecX);

rotMat = [newVecX; newVecY; newVecZ];

end

