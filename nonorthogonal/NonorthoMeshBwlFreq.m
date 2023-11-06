function bwl = NonorthoMeshBwlFreq(a1, a2, na1, na2, n1, n2, prop)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

a3 = [0, 0, 1];
[b1, b2, ~] = DirectToReciprocal(a1, a2, a3);
b1 = b1 * n1 / na1 / (2 * pi);
b2 = b2 * n2 / na2 / (2 * pi);

b1b2 = norm(cross(b1, b2));
projB1 = b1b2 / norm(b2);
projB2 = b1b2 / norm(b1);

bwl = 1/2 * min(projB1, projB2) * prop;

end